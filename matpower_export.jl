# --------------------------------------------------------------------------- #
# EXPORT -- write the reduced network back out as a MATPOWER .m file, so it
# can be loaded directly by MATPOWER or PowerModels, independent of this repo.
#
# The reduction here is pure bus contraction: an internal line's two endpoints
# become electrically identical, so nothing about a SURVIVING line's own data
# ever changes. That is what makes this safe to build on top of
# PowerModels.export_matpower rather than writing MATPOWER syntax by hand --
# the only edits below are which bus_i a load/shunt/generator/branch endpoint
# points at. Resistance, susceptance, tap ratio, phase shift, rating and
# status of every surviving branch are untouched, because the branch dict
# entry itself is never modified, only re-pointed. export_matpower then does
# the per-unit -> physical-unit conversion and the load/shunt aggregation by
# bus itself -- the same code path any ordinary MATPOWER file already goes
# through, not something written new for this purpose.
#
# CAVEAT worth carrying into any README next to an exported file: this whole
# codebase's own DC model (tnr_preprocessing.jl's build_tx_case) already
# ignores tap ratio and phase shift when computing flows -- Dx = 1/br_x only.
# That is a pre-existing simplification of the whole pipeline, not something
# this export introduces. A tap/shift-aware DC or AC solve elsewhere may
# therefore see slightly different flows than what this repo validated.
# --------------------------------------------------------------------------- #

using PowerModels

"""
    export_reduced_matpower(casefile, c, r, outfile) -> outfile

Write the network reduced by clustering `r` (full-network granularity: an
r_display from a Kron run, never r_boundary -- see kron_unfold.jl) as a
MATPOWER .m file at `outfile`.

Retained buses keep their ORIGINAL external bus numbers (traceable against
bus_mapping.csv); eliminated buses are simply absent. Every load and shunt
that belonged to an eliminated bus is redirected onto its cluster's
representative -- not summed by hand here, but by export_matpower itself once
`load_bus`/`shunt_bus` point at the representative, exactly mirroring how
`_dcopf_on_partition` (tnr_postprocessing.jl) already aggregates load onto
representative buses for the validated reduced-network benchmark. Generators
are never aggregated, only `gen_bus` is redirected -- each keeps its own
limits and cost.

Errors on storage/switch/dcline components: none of the bundled cases carry
them, and silently dropping such a component's bus reference on a case that
DOES have one would misrepresent the network rather than merely shrink it.
"""
function export_reduced_matpower(casefile::AbstractString, c, r, outfile::AbstractString)
    raw = PowerModels.parse_file(casefile)
    for comp in ("storage", "switch", "dcline")
        isempty(raw[comp]) ||
            error("export_reduced_matpower: $casefile has $(length(raw[comp])) " *
                  "'$comp' entries, which this exporter does not redirect on " *
                  "bus contraction -- extend it before using it on this case.")
    end

    base = c.base
    # c.bus_ids is the COMPACT number PowerModels.make_basic_network assigns
    # internally (1..N, gaps removed) -- not the original file's bus_i
    # whenever the file itself has gaps or non-contiguous numbering (e.g.
    # pglib_opf_case300's 9000-series buses). make_basic_network keeps the
    # true original number in each bus's source_id, so recover it from a
    # matching basic-network parse rather than trusting bus_ids to already be
    # it. For an already-contiguous file (case118, ACTIVSg*) this is a no-op:
    # compact number == original number.
    basic = make_basic_network(PowerModels.parse_file(casefile))
    compact_to_true = Dict{Int,Int}(b["bus_i"] => b["source_id"][2] for (_, b) in basic["bus"])
    true_id = [compact_to_true[c.bus_ids[i]] for i in 1:base.N]
    rep_id = Dict{Int,Int}(true_id[i] => true_id[r.rep_of[i]] for i in 1:base.N)

    for k in collect(keys(raw["bus"]))
        b = raw["bus"][k]
        rep_id[b["bus_i"]] == b["bus_i"] || delete!(raw["bus"], k)
    end
    for (_, ld) in raw["load"]
        ld["load_bus"] = rep_id[ld["load_bus"]]
    end
    for (_, sh) in raw["shunt"]
        sh["shunt_bus"] = rep_id[sh["shunt_bus"]]
    end
    for (_, g) in raw["gen"]
        g["gen_bus"] = rep_id[g["gen_bus"]]
    end
    # A branch survives iff its endpoints land in different clusters after
    # remapping -- true regardless of the branch's own original in-service
    # status, and correctly drops an already-out-of-service branch whose
    # endpoints happen to merge through some OTHER line, avoiding a stray
    # f_bus == t_bus row.
    for k in collect(keys(raw["branch"]))
        br = raw["branch"][k]
        fr, to = rep_id[br["f_bus"]], rep_id[br["t_bus"]]
        if fr == to
            delete!(raw["branch"], k)
        else
            br["f_bus"], br["t_bus"] = fr, to
        end
    end

    mkpath(dirname(outfile))
    PowerModels.export_matpower(outfile, raw)
    return outfile
end
