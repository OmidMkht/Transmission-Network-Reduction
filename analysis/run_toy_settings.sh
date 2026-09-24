#!/bin/bash
# --------------------------------------------------------------------------- #
# Settings sweep for the three reduction approaches on a toy case.
#
#   bash analysis/run_toy_settings.sh 14bus
#   bash analysis/run_toy_settings.sh 6busww
#
# Window: a three-month quarter (Sep-Oct-Nov, 2184 h). The pipeline's seed
# selection ranks every hour of that quarter by total |line flow| and takes the
# 20 highest as TRAINING; the other ~2164 are the test set. So the training set
# is the high-flow hours of the whole window, not of one month.
#
#   KKT     cost cap 1 / 2 / 3 %
#   Greedy  cost cap 1 / 2 / 3 % (matched to KKT), plus a loose flow tolerance
#   Proxy   flow window eps 0.10 / 0.25 / 0.50, plus two near_limit settings.
#           eps below 0.10 produced nothing at all on either case -- the proxy
#           only ever returned the radial merge its preprocessing finds free.
#
# No hop cap, no line budget, no cluster-size cap anywhere: the toys are small
# enough to solve outright. The proxy also runs with lmp_separation=false and
# cycle_cuts=() so no heuristic cut restricts any approach; kkt already passes
# cycle_cut_lens=() and defaults lmp_separation to false, and greedy has none.
#
# output_dir is set explicitly for kkt and greedy: their run_tag does NOT
# include cost_cap, so every setting would otherwise overwrite the same folder.
# --------------------------------------------------------------------------- #
set -uo pipefail
cd "$(dirname "$0")/.."

CASE="${1:-14bus}"
case "$CASE" in
  6busww) ID=case6ww  ;;
  14bus)  ID=case14toy ;;
  *) echo "unknown case $CASE (6busww | 14bus)" >&2; exit 2 ;;
esac

MONTHS="[9,10,11]"
W="case=$ID demands=hourly month=$MONTHS horizon_start_day=1 horizon_days=92 n_demands=20"
TL=300
J="julia --project=. --startup-file=no"

# Every run is wall-clocked and the elapsed seconds written next to its output,
# so the review can put a runtime column beside the quality columns. Wall time
# is the honest operational number -- it is what "time to find the reduction"
# means -- and it is the only measure available for all three, since the proxy
# records no solve time of its own. It includes ~20 s of Julia start-up.
timed() {            # timed <outdir> <cmd...>
    local out="$1"; shift
    local t0=$SECONDS
    "$@" 2>&1 | grep -E "^kkt |^greedy  |^proxy  " | tail -1
    mkdir -p "$out"
    echo $(( SECONDS - t0 )) > "$out/runtime_s.txt"
}

# start from a clean sheet so the review table shows this sweep only
rm -rf "outputs/kkt/$ID" "outputs/greedy/$ID" "outputs/proxy/$ID"

for cap in 1 2 3; do
    echo "=== KKT  cost_cap=${cap}% ==="
    timed "outputs/kkt/$ID/cap${cap}pct" $J kkt/run_kkt.jl $W cost_cap=$cap time_limit=$TL "output_dir=outputs/kkt/$ID/cap${cap}pct"
done

for cap in 1 2 3; do
    echo "=== Greedy  cost_cap=${cap}% ==="
    timed "outputs/greedy/$ID/cap${cap}pct" $J greedy/run_greedy.jl $W cost_cap=$cap time_limit=$TL "output_dir=outputs/greedy/$ID/cap${cap}pct"
done

echo "=== Greedy  cost_cap=3%, loose flow_tol=1e-4 ==="
timed "outputs/greedy/$ID/cap3pct_tol1e-4" $J greedy/run_greedy.jl $W cost_cap=3 flow_tol=1e-4 time_limit=$TL "output_dir=outputs/greedy/$ID/cap3pct_tol1e-4"

# eps alone, at the default near_limit = 0.8
for eps in 0.10 0.25 0.50; do
    tag="eps${eps//./p}"
    echo "=== Proxy  eps=$eps  near_limit=0.8 ==="
    timed "outputs/proxy/$ID/$tag" $J proxy/run_proxy.jl $W eps=$eps time_limit=$TL plots=false lmp_separation=false "cycle_cuts=()" "output_dir=outputs/proxy/$ID/$tag"
done

# near_limit decides how many congested lines are held out of every cluster,
# and on case300 it -- not eps -- turned out to be the binding knob.
for nl in 0.95 0.99; do
    tag="eps0p25_nl${nl//./p}"
    echo "=== Proxy  eps=0.25  near_limit=$nl ==="
    timed "outputs/proxy/$ID/$tag" $J proxy/run_proxy.jl $W eps=0.25 near_limit=$nl time_limit=$TL plots=false lmp_separation=false "cycle_cuts=()" "output_dir=outputs/proxy/$ID/$tag"
done

echo
echo "=== review ==="
$J analysis/review_toy.jl case=$CASE month=9,10,11 n_train=20 all_settings=true 2>&1 | grep -vE "LicenseID|Academic|Username|^Set parameter|redispatch|^Automatic"
