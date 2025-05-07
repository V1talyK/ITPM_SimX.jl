function save_rgm2file(tag, rsl, grd, prp, gdm_prop, wxy)
    rt = joinpath(Base.source_dir(),"v4")
    tlb = Vector(undef,0)
    for iw = 1:nw
        for t = 1:nt
            qp = ifelse(rsl.qw[iw,t] >= 0.0, rsl.qw[iw,t], 0.0)
            qi = ifelse(rsl.qw[iw,t] <= 0.0, rsl.qw[iw,t], 0.0)
            qo = ifelse(rsl.qw[iw,t] >= 0.0, rsl.qw[iw,t], 0.0)
            opra_i = qi > 0.0
            opra_p = qp > 0.0
            tmp = Vector{Any}(undef, 0)
            append!(tmp, iw, t, gdm_prop.dt*t, qp, qo, qp.-qo, qi, rsl.pw[iw,t], rsl.ppl[iw,t], opra_p, opra_i)
            push!(tlb, tmp)
        end
    end


    head = ["id","numb","date","liquid","oil","water","injection","bhp","cellp_m","opra_p","opra_i"]
    pth = joinpath(rt,"script_$(tag)_2500.csv")
        write_to_csv(pth, tlb, head)
    
    prm = Dict("bnd"=>[grd.X[grd.λbi.-(grd.nc-grd.nc0)], grd.Y[grd.λbi.-(grd.nc-grd.nc0)]],
               "wxy"=>wxy,
               "wkp"=>prp.kp[getindex.(well,1)],
               "whe"=>prp.he[getindex.(well,1)],
               "kp"=>round.(prp.kp, sigdigits=3),
               "he"=>round.(prp.he, sigdigits=3),
               "mp"=>round.(prp.mp, sigdigits=3),
               "prop"=>gdm_prop,
               "2f"=>Dict("mu_w"=>1, "mu_o"=>3, "n_w"=>2, "n_o"=>2))
    
    pth = joinpath(rt ,"script_$(tag)_2500_2f.json")
        write_to_json(pth, prm)
    
        return nothing
end


function make_dsp(nw, rsl)
    dsp = zeros(nw)
    for iw = 1:nw
        qwi = rsl.qw[iw, :]
        mima = maximum(qwi) .- minimum(qwi)
        if mima == 0

        else
            qwi = (qwi .- minimum(qwi)) ./ mima
        end
        tmp = abs.(diff(qwi))
        mtmp = median(tmp)
        ia = findall(tmp .> mtmp)
        ia = vcat(1, ia .+ 1)
        qr = zeros(0)
        for i = 1:length(ia)
            if i == length(ia)
                li = nt
            else
                li = ia[i+1]
            end
            push!(qr, mean(qwi[ia[i]:li]))
        end
        # dsp[iw] = sqrt.(mean((mean(qr).-qr).^2))
        # gap = cut_by_production(qwi.+1)
        # gap[end] = gap[end]-1
        # qr = zeros(0)
        # for i = 1:length(gap)-1
        #     push!(qr, mean(qwi[gap[i]:gap[i+1]]))
        # end
        dsp[iw] = sqrt.(mean((mean(qr) .- qr) .^ 2))
    end
    return dsp
end

function fix_qw(qw, uf, pw)
  flag = true
  k = 0
  kmax = 20;
  while flag
    k += 1
    rsl = sim_calc(qw=qw, uf=uf, pw=pw)
    ia = findall(rsl.pw .< 0.05)
    ib = findall(rsl.pw .> 25)
    flag1 = (&)(k < kmax, length(ia) > 0)
    if flag1
      tmp = 0.05.-rsl.pw[ia];
      max_er_pw = maximum(tmp);
      alp = tmp./max_er_pw;
      qw[ia] .*= 1.0 .- 0.3.*alp
    end
    flag2 = (&)(k < kmax, length(ib) > 0)
    if flag2
      qw[ib] .*= 0.9
    end
    println(length(ia))
    flag = flag1 | flag2
  end
  return qw
end

sim_calc, cIWC = make_sim(grd, gdm_prop, well, prp, nt);

tag = "s0.1"
qw = ones(Float32, nw, nt);
#qw, _uf = SimScriptTool.gen_real_rand_qw(nw, nt; mult = 1.0)
qw[iw_inj,:] .*= -1
qw = SimScriptTool.scaleRateByGeoProp(qw, prp.Vp, prp.kp, prp.he, gdm_prop.dt, gdm_prop.P0, grd.ds)

pw = 2*ones(nw, nt);
uf =  falses(nw, nt);
uf[iw_inj,:] .= true
pw[iw_inj,:] .= 12
qw = fix_qw(qw, uf, pw);
qw .= mean(qw, dims = 2);
rsl = sim_calc(qw = qw, uf = uf, pw = pw)
save_rgm2file(tag, rsl, grd, prp, gdm_prop, wxy)
sum(make_dsp(nw, rsl))

tag = "s0.2"
for (k,v) in enumerate(zip(Iterators.partition(169:nt,floor(Int64, (nt-168)/nw)), Iterators.cycle(1:nw)))
    #println(k," ",v)
    qw[v[2],v[1]] .*= 0.5;
end
rsl = sim_calc(qw = qw, uf = uf, pw = pw)
save_rgm2file(tag, rsl, grd, prp, gdm_prop, wxy)
sum(make_dsp(nw, rsl))

tag = "s0.3"
qw = ones(Float32, nw, nt);
#qw, _uf = SimScriptTool.gen_real_rand_qw(nw, nt; mult = 1.0)
#qw[iw_inj,:] .*= -1
qw = SimScriptTool.scaleRateByGeoProp(qw, prp.Vp, prp.kp, prp.he, gdm_prop.dt, gdm_prop.P0, grd.ds)

pw = 2*ones(nw, nt);
uf =  falses(nw, nt);
qw = fix_qw(qw, uf, pw);
qw .= mean(qw, dims = 2).*0.95
rsl = sim_calc(qw = qw, uf = uf, pw = pw)
save_rgm2file(tag, rsl, grd, prp, gdm_prop, wxy)
sum(make_dsp(nw, rsl))


tag = "s0.4"
for (k,v) in enumerate(zip(Iterators.partition(169:nt,floor(Int64, (nt-168)/nw)), Iterators.cycle(1:nw)))
  #println(k," ",v)
  qw[v[2],v[1]] .*= 0.5;
end
rsl = sim_calc(qw = qw, uf = uf, pw = pw)
save_rgm2file(tag, rsl, grd, prp, gdm_prop, wxy)
sum(make_dsp(nw, rsl))


tag = "s1.1"
qw = ones(Float32, nw, nt);
qw, _uf = SimScriptTool.gen_real_rand_qw(nw, nt; mult = 1.0)
qw[iw_inj,:] .*= -1
qw = SimScriptTool.scaleRateByGeoProp(qw, prp.Vp, prp.kp, prp.he, gdm_prop.dt, gdm_prop.P0, grd.ds)
pw = 2*ones(nw, nt);
uf =  falses(nw, nt);
uf[iw_inj,:] .= true
pw[iw_inj,:] .= 14
qw = fix_qw(qw, uf, pw)
rsl = sim_calc(qw = qw, uf = uf, pw = pw)
save_rgm2file(tag, rsl, grd, prp, gdm_prop, wxy)
sum(make_dsp(nw, rsl))


tag = "s1.2"
qw = SimScriptTool.gen_rand_qw(nw, nt)
qw[iw_inj,:] .*= -1
qw = SimScriptTool.scaleRateByGeoProp(qw, prp.Vp, prp.kp, prp.he, gdm_prop.dt, gdm_prop.P0, grd.ds)

pw = 2*ones(nw, nt);
uf =  falses(nw, nt);
uf[iw_inj,:] .= true
pw[iw_inj,:] .= 14
qw = fix_qw(qw, uf, pw)
rsl = sim_calc(qw = qw, uf = uf, pw = pw)
save_rgm2file(tag, rsl, grd, prp, gdm_prop, wxy)
sum(make_dsp(nw, rsl))

tag = "s1.3"
qw = SimScriptTool.gen_periodic(nw, nt, cnt = 5)
qw[iw_inj,:] .*= -1
qw = SimScriptTool.scaleRateByGeoProp(qw, prp.Vp, prp.kp, prp.he, gdm_prop.dt, gdm_prop.P0, grd.ds)

pw = 2*ones(nw, nt);
uf =  falses(nw, nt);
uf[iw_inj,:] .= true
pw[iw_inj,:] .= 14
qw = fix_qw(qw, uf, pw)
rsl = sim_calc(qw = qw, uf = uf, pw = pw)
save_rgm2file(tag, rsl, grd, prp, gdm_prop, wxy)
sum(make_dsp(nw, rsl))


tag = "s1.4"
fl = SimScriptTool.gen_rgm(nw);
qw = ones(Float32, nw, nt);
qw[iw_inj,:] .*= -1
qw[.!fl[:,1:nt]].*=0.1;
qw = SimScriptTool.scaleRateByGeoProp(qw, prp.Vp, prp.kp, prp.he, gdm_prop.dt, gdm_prop.P0, grd.ds)

pw = 2*ones(nw, nt);
uf =  falses(nw, nt);
uf[iw_inj,:] .= true
pw[iw_inj,:] .= 14
qw = fix_qw(qw, uf, pw)
rsl = sim_calc(qw = qw, uf = uf, pw = pw)
save_rgm2file(tag, rsl, grd, prp, gdm_prop, wxy)
sum(make_dsp(nw, rsl))


tag = "s2.1"
qw = ones(Float32, nw, nt);
#qw, _uf = SimScriptTool.gen_real_rand_qw(nw, nt; mult = 1.0)
qw[iw_inj,:] .*= -1
qw = SimScriptTool.scaleRateByGeoProp(qw, prp.Vp, prp.kp, prp.he, gdm_prop.dt, gdm_prop.P0, grd.ds)

pw = 2*ones(nw, nt);
uf =  falses(nw, nt);
uf[iw_inj,:] .= true
pw[iw_inj,:] .= 12
qw = fix_qw(qw, uf, pw)
qw .= mean(qw, dims = 2)
for (k, v) in enumerate(zip(Iterators.partition(1:nt, ceil(Int64, 168/nw)), Iterators.cycle(1:nw)))
  #enumerate(Iterators.partition(1:nt,ceil(Int64,nt/nw)))
  #println(v)
  qw[v[2],v[1]] .= 0f0; 
  uf[v[2],v[1]] .= false;
end
qw = fix_qw(qw, uf, pw)
rsl = sim_calc(qw = qw, uf = uf, pw = pw)

save_rgm2file(tag, rsl, grd, prp, gdm_prop, wxy)
sum(make_dsp(nw, rsl))


tag = "s2.2"
qw = ones(Float32, nw, nt);
#qw, _uf = SimScriptTool.gen_real_rand_qw(nw, nt; mult = 1.0)
qw[iw_inj,:] .*= -1
qw = SimScriptTool.scaleRateByGeoProp(qw, prp.Vp, prp.kp, prp.he, gdm_prop.dt, gdm_prop.P0, grd.ds)

pw = 2*ones(nw, nt);
uf =  falses(nw, nt);
uf[iw_inj,:] .= true
pw[iw_inj,:] .= 12
qw = fix_qw(qw, uf, pw)
qw .= mean(qw, dims = 2)
for (k, v) in enumerate(zip(Iterators.partition(1:nt, ceil(Int64, 168/nw)), Iterators.cycle(1:nw)))
  qw[v[2],v[1]] .*= 0.5f0; 
  qw[v[2],v[1]] .= false;
end
rsl = sim_calc(qw = qw, uf = uf, pw = pw)

save_rgm2file(tag, rsl, grd, prp, gdm_prop, wxy)
sum(make_dsp(nw, rsl))