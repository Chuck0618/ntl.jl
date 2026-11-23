using NTL, Statistics
Base.@noinline function checksum(v)
    s=0
    for x in v
        s+=Int(rep(x))
    end
    s
end
function report()
    F=PrimeField(97)
    input=mod.(0:99999,97)
    packed=packedfieldvec(F,input)
    plain=F.(input)
    @assert checksum(packed)==checksum(plain)==sum(input)
    for (name,v,payload) in (("packed",packed,sizeof(packed.data)),("ordinary",plain,sizeof(plain)))
        times=Float64[]
        for _ in 1:51
            start=time_ns()
            observed=checksum(Base.inferencebarrier(v))
            elapsed=time_ns()-start
            @assert observed==sum(input)
            push!(times,elapsed/1e9)
        end
        println(name," payload_bytes=",payload," median_read_us=",round(median(times)*1e6;digits=2))
    end
end
report()
