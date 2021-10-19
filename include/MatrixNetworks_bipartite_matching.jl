
function bipartite_matching_primal_dual(rp::Vector{Int64}, ci::Vector{Int64}, 
                    ai::Vector{T}, m::Int64, n::Int64) where T
    
    # variables used for the primal-dual algorithm
    # normalize ai values # updated on 2-19-2019
    #ai ./= maximum(abs.(ai))
    alpha=zeros(Float64,m)
    bt=zeros(Float64,m+n)#beta
    queue=zeros(Int64,m)
    t=zeros(Int64,m+n)
    match1=zeros(Int64,m)
    match2=zeros(Int64,m+n)
    tmod = zeros(Int64,m+n)
    ntmod=0
    
    # initialize the primal and dual variables
    for i=1:m
        for rpi=rp[i]:rp[i+1]-1
            if ai[rpi] > alpha[i]
               alpha[i]=ai[rpi]
            end
        end
    end
    
    # dual variables (bt) are initialized to 0 already
    # match1 and match2 are both 0, which indicates no matches
    
    i=1
    while i<=m
        for j=1:ntmod
            t[tmod[j]]=0
        end
        ntmod=0
        # add i to the stack
        head=1
        tail=1
        queue[head]=i
        while head <= tail && match1[i]==0
            k=queue[head]
            for rpi=rp[k]:rp[k+1]-1
                j = ci[rpi]
                if ai[rpi] < alpha[k] + bt[j] - 1e-8
                    continue
                end # skip if tight
                if t[j]==0
                    tail=tail+1
                    if tail <= m
                        queue[tail]=match2[j]
                    end
                    t[j]=k
                    ntmod=ntmod+1
                    tmod[ntmod]=j
                    if match2[j]<1
                        while j>0
                            match2[j]=t[j]
                            k=t[j]
                            temp=match1[k]
                            match1[k]=j
                            j=temp
                        end
                        break
                    end
                end
            end
            head=head+1
        end
        if match1[i] < 1
            theta=Inf
            for j=1:head-1
                t1=queue[j]
                for rpi=rp[t1]:rp[t1+1]-1
                    t2=ci[rpi]
                    if t[t2] == 0 && alpha[t1] + bt[t2] - ai[rpi] < theta
                        theta = alpha[t1] + bt[t2] - ai[rpi]
                    end
                end
            end
            for j=1:head-1
                alpha[queue[j]] -= theta
            end
            for j=1:ntmod
                bt[tmod[j]] += theta
            end
            continue
        end
        i=i+1
    end
    val=0
    for i=1:m
        for rpi=rp[i]:rp[i+1]-1
            if ci[rpi]==match1[i]
                val=val+ai[rpi]
            end
        end
    end
    noute = 0
    for i=1:m
        if match1[i]<=n
            noute=noute+1
        end
    end

    M_output = MatrixNetworks.Matching_output(m,n,val,noute,match1)
    return M_output
end
