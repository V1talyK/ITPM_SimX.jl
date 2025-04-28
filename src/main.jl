push!(LOAD_PATH, Base.source_dir())

using ITPM_SimX
using UnicodePlots, StatsBase
using Revise

@info "Тестирование прямого расчёта симулятора"
grd, gdm_prop, prp, xy, nt = make_gdm(kp_init = 0.5, 
                                    λb = 0.0)
grd, prp = aq_extend(grd, gdm_prop, prp; prm=Dict("nor" => 1,
                                                  "lat" => 1,
                                                  "vol" => 100))

wxy9 = collect(zip(1:9,collect(Iterators.product(xy[1],xy[2]))[:]))
insert!(wxy9, 6, (5, (500, 600)))
well = make_well(wxy9,grd)
nw = length(unique(getindex.(well,2)))

sim_calc, cIWC = make_sim(grd, gdm_prop, well, prp, nt)

qw = rand(-1:0.1:1, nw, nt);
qw[[1,3,7,9],:] .= -abs.(qw[[1,3,7,9],:]).-1;
qw[[2,4,5,6,8],:] .= abs.(qw[[2,4,5,6,8],:].+1)

rsl = sim_calc(qw = qw)
wc0 = ones(nw, nt)
wc0[1, :] .= 2
rsl1 = sim_calc(qw = qw, wc = wc0)
@btime sim_calc(qw = $qw)
@profiler sim_calc(qw = qw)

lineplot(rsl.ppl[1,:])|>println
heatmap(reshape(rsl.PM[1:grd.nc0,1], floor(Int64, sqrt(grd.nc0)), floor(Int64,sqrt(grd.nc0))))|>println
lineplot(mean(rsl.PM[grd.nc0+1:end,:], dims=1)[:])|>println

uf = falses(nw, nt);
uf[[2,4,5,6,8],:] .= true;
qw[[1,3,7,9],:] .= -2;
pw = zeros(nw, nt);
pw[[2,4,5,6,8],:] .= 5f0;
#pw[[2,4,5,6,8],:] .= rsl.pw[[2,4,5,6,8],:];

rsl = sim_calc(qw = qw, uf = uf, pw = pw)
