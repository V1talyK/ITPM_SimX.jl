using Plots, SimScriptTool
using JLD2, FileIO
rsrc = Base.source_dir()
DD = load(joinpath(rsrc,"grd.jld2"))
DD["XY"] = DD["XY"] .- minimum(DD["XY"],dims=1)

grd, gdm_prop, prp, nt = make_gdmV(;he_init = 40.,
                                     kp_init = 100,
                                     mp_init = 0.2,
                                     nt_init = 120,
                                     DD = DD,
                                     Paq = 8.6,
                                     bet = 1e-4,
                                     Paq = 8,
                                     λb = 2.0)


wxy101 = [DD["XY"][DD["Won"][iw.==DD["Won"][:,2],1], :] for iw in DD["wi"]]
well = collect(zip(DD["Won"][:,1], indexin(DD["Won"][:,2], DD["wi"])))
nw = length(unique(getindex.(well,2)))

sim_calc, cIWC = make_sim(grd,gdm_prop, well, prp, nt);
qw = rand(-1:0.1:1, nw, nt);
qw, _uf = SimScriptTool.gen_real_rand_qw(nw, nt; mult = 1.0)
qw .= abs.(qw)
qw[iw_inj,:] .*= -1
qw = SimScriptTool.scaleRateByGeoProp(qw, prp.Vp, prp.kp, prp.he, gdm_prop.dt, gdm_prop.P0, grd.ds)

pw = 2*ones(nw, nt);
uf =  falses(nw, nt)


rsl = sim_calc(qw = qw, uf = uf, pw = pw)
sum(rsl.qw[rsl.qw.>0])*30.4./sum(prp.Vp)
iw = 4
  plt = plot(rsl.ppl[iw,:])
  plot!(plt, rsl.pw[iw,:])

plot(rsl.qw[iw,:])
qp = sum.(filter.(x->x>0, eachcol(rsl.qw)))
qi = sum.(filter.(x->x<0, eachcol(rsl.qw)))
plt = plot(qp)
  plot(plt, -qi)
using StatsBase
plot(mean(rsl.ppl, dims =1)[:])

plt = plot()
  plt2 = plot()
for (k,v) in enumerate([g1, g2, g3, g4])
  qp = sum.(filter.(x->x>0, eachcol(rsl.qw[v,:])))
  qi = sum.(filter.(x->x<0, eachcol(rsl.qw[[iw_inj[k]],:])))
  pplG = mean.(filter.(x->x>0, eachcol(rsl.ppl[v,:])))
  plot!(plt2, pplG)
  plot!(plt, qp)
  plot!(plt, -qi)
end
display(plt)
display(plt2)

    #
using Dates

pth = "/home/lik/rgm_v1"
tlb = []
vd = range(Date("2000-01-01"), step = Dates.Month(1), length=nt)
for iw = 1:nw
  for t = 1:nt
    prod = rsl.qw[iw,t]>0 ? rsl.qw[iw,t] : 0.0
    inj = rsl.qw[iw,t]<0 ? -rsl.qw[iw,t] : 0.0
    push!(tlb,[iw, vd[t], rsl.pw[iw,t]*10,prod, inj])
  end
end
write_to_csv(pth, tlb, ["well", "date", "pw, атм.", "prod, м3/сут.", "inj, м3/сут."])
