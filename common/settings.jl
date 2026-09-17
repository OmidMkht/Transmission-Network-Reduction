# Settings block + command-line overrides for the runners.
#
#   julia --project=. greedy/run_greedy.jl case=case300 hop_cap=5 time_limit=1800
#
# Values: numbers, true/false, nothing, :symbol or a bare word, [vectors], (tuples).
# A bare word becomes a Symbol when the default is a Symbol, otherwise a String.

module Settings

export apply_overrides, print_settings, save_settings

const NOT_LITERAL = gensym(:not_literal)

function literal(ex)
    ex isa Number && return ex
    ex isa AbstractString && return String(ex)
    ex isa QuoteNode && ex.value isa Symbol && return ex.value
    ex === :nothing && return nothing
    ex === :true && return true
    ex === :false && return false
    if ex isa Expr && ex.head in (:vect, :tuple)
        vals = map(literal, ex.args)
        any(x -> x === NOT_LITERAL, vals) && return NOT_LITERAL
        return ex.head === :vect ? map(identity, vals) : Tuple(vals)
    end
    return NOT_LITERAL
end

function parse_value(raw::AbstractString, default)
    s = strip(raw)
    s in ("nothing", "none") && return nothing
    ex = try
        Meta.parse(s)
    catch
        NOT_LITERAL
    end
    v = literal(ex)
    v === NOT_LITERAL || return v
    default isa Symbol && return Symbol(lstrip(s, ':'))
    return String(s)
end

"Merge `key=value` arguments into the settings, rejecting unknown keys."
function apply_overrides(settings::NamedTuple, args=ARGS)
    isempty(args) && return settings
    updates = Pair{Symbol,Any}[]
    for a in args
        occursin('=', a) || error("expected key=value, got \"$a\"")
        k, v = split(a, '='; limit=2)
        key = Symbol(strip(k))
        haskey(settings, key) || error("unknown setting \"$key\". Valid: " *
                                       join(string.(keys(settings)), ", "))
        push!(updates, key => parse_value(v, settings[key]))
    end
    return merge(settings, NamedTuple(updates))
end

function print_settings(settings::NamedTuple; io::IO=stdout)
    w = maximum(length(string(k)) for k in keys(settings))
    for (k, v) in pairs(settings)
        println(io, "  ", rpad(string(k), w), " = ", repr(v))
    end
end

function save_settings(path::AbstractString, settings::NamedTuple)
    mkpath(dirname(path))
    open(path, "w") do io
        print_settings(settings; io)
    end
end

end # module
