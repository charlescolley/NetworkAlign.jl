using MatrixNetworks
include("normout.jl")
"""
ISORANK
-------
    solve an overlap matching problem with IsoRank
"""

function isorank{T}(S::SparseMatrixCSC{T,Int64},w::Vector{Float64},
                 a::Int64,b::Int64,li::Vector{Int64},lj::Vector{Int64},
                 alpha::Float64,rtype::Int64,tol::Float64,
                 maxit::Int64,verbose::Bool,P::SparseMatrixCSC{Float64,Int64})
  csum = sum_kbn(w)
  v = w./csum
  assert(all(v.>=0))
  n = size(P,1)

  allstats = !isempty(li) && !isempty(lj)
  if allstats
    rhistsize = 5
    M_setup = MatrixNetworks.bipartite_matching_setup(w,li,lj,maximum(li),maximum(lj))
    tripi = M_setup.tripi
    matm = M_setup.m
    matn = M_setup.n
    rp = M_setup.rp
    ci = M_setup.ci
    mperm = tripi[tripi.>0] # a permutation for the matching problem
  else
    rhistsize = 1
  end
  r = alpha
  x = zeros(Float64,n) + v
  delta = 2
  iter = 0
  reshist = zeros(Float64,maxit,rhistsize)
  xbest = x
  fbest = 0
  fbestiter = 0
  if verbose && allstats #print the header
    @printf("%5s   %4s   %8s   %7s %7s %7s %7s\n","best","iter","pr-delta","obj","weight","card","overlap")
  elseif verbose
    @printf("%4s   %8s\n","iter","delta")
  end
  while iter < maxit && delta > tol
    y = r * (P' * x)
    omega = sum_kbn(x) - sum_kbn(y)
    y = y + omega * v
    delta = vecnorm(x-y,1)
    reshist[iter+1] = delta
    iter = iter + 1
    x = y * (1/sum_kbn(y))
    if allstats
      if rtype == 1
        xf = x
      elseif rtype == 2
        xf = a*v + b/2*(S*x) #add v to preserve scale
      end
      ai = zeros(Float64,length(tripi))
      ai[tripi.>0] = xf[mperm]
      M_output = MatrixNetworks.bipartite_matching_primal_dual(rp,ci,ai,matm,matn)
      ma = M_output.cardinality
      #get match:
      mi_int = MatrixNetworks.edge_indicator(M_output,li,lj)
      val = dot(w,mi_int)
      overlap = dot(mi_int,(S*mi_int)/2)
      f = a*val + b*overlap
      if f > fbest
        xbest = x
        fbest = f
        fbestiter = iter
        itermark = "*"
      else
        itermark = " "
      end
      if verbose && allstats
        @printf("%5s   %4i   %8.1e   %5.2f %7.2f %7i %7i\n", itermark, iter, delta, f, val, ma, overlap)
        reshist[iter,2:end] = [a*val + b*overlap, val, ma, overlap]
      elseif verbose
        @printf("%4i    %8.1e\n", iter, delta)
      end
    end
  end
  flag = delta>tol
  reshist = reshist[1:iter,:]
  if allstats
    x=xbest
  end
  return (x,flag,reshist)
end

function isorank{T}(S::SparseMatrixCSC{T,Int64},w::Vector{Float64},
                 a::Int64,b::Int64,li::Vector{Int64},lj::Vector{Int64},
                 alpha::Float64,rtype::Int64,tol::Float64,
                 maxit::Int64,verbose::Bool)
  ss = sum(S,2)
  P = sparse(S./ss)
  return isorank(S,w,a,b,li,lj,alpha,rtype,tol,maxit,verbose,P)
end

function isorank{T}(S::SparseMatrixCSC{T,Int64},w::Vector{Float64},
                 a::Int64,b::Int64,li::Vector{Int64},lj::Vector{Int64},
                 alpha::Float64,rtype::Int64,tol::Float64,
                 maxit::Int64)
  verbose = true
  ss = sum(S,2)
  P = sparse(S./ss)
  return isorank(S,w,a,b,li,lj,alpha,rtype,tol,maxit,verbose,P)
end

function isorank{T}(S::SparseMatrixCSC{T,Int64},w::Vector{Float64},
                 a::Int64,b::Int64,li::Vector{Int64},lj::Vector{Int64},
                 alpha::Float64,rtype::Int64,tol::Float64)
  maxit = 100
  verbose = true
  ss = sum(S,2)
  P = sparse(S./ss)
  return isorank(S,w,a,b,li,lj,alpha,rtype,tol,maxit,verbose,P)
end

function isorank{T}(S::SparseMatrixCSC{T,Int64},w::Vector{Float64},
                 a::Int64,b::Int64,li::Vector{Int64},lj::Vector{Int64},
                 alpha::Float64,rtype::Int64)
  tol = 1e-12
  maxit = 100
  verbose = true
  ss = sum(S,2)
  P = sparse(S./ss)
  return isorank(S,w,a,b,li,lj,alpha,rtype,tol,maxit,verbose,P)
end

function isorank{T}(S::SparseMatrixCSC{T,Int64},w::Vector{Float64},
                 a::Int64,b::Int64,li::Vector{Int64},lj::Vector{Int64},
                 alpha::Float64)
  rtype = 2
  tol = 1e-12
  maxit = 100
  verbose = true
  ss = sum(S,2)
  P = sparse(S./ss)
  return isorank(S,w,a,b,li,lj,alpha,rtype,tol,maxit,verbose,P)
end

function isorank{T}(S::SparseMatrixCSC{T,Int64},w::Vector{Float64},
                 a::Int64,b::Int64,li::Vector{Int64},lj::Vector{Int64})
  alpha = b/(a+b)
  rtype = 2
  tol = 1e-12
  maxit = 100
  verbose = true
  P = normout(S)
  return isorank(S,w,a,b,li,lj,alpha,rtype,tol,maxit,verbose,P)
end

function isorank{T}(S::SparseMatrixCSC{T,Int64},w::Vector{Float64},
                 a::Int64,b::Int64)
  li = []
  lj = []
  alpha = b/(a+b)
  rtype = 2
  tol = 1e-12
  maxit = 100
  verbose = true
  ss = sum(S,2)
  P = sparse(S./ss)
  return isorank(S,w,a,b,li,lj,alpha,rtype,tol,maxit,verbose,P)
end

function isorank{T}(S::SparseMatrixCSC{T,Int64},w::Vector{Float64},
                 a::Int64)
  b = 1
  li = []
  lj = []
  alpha = b/(a+b)
  rtype = 2
  tol = 1e-12
  maxit = 100
  verbose = true
  ss = sum(S,2)
  P = sparse(S./ss)
  return isorank(S,w,a,b,li,lj,alpha,rtype,tol,maxit,verbose,P)
end

function isorank{T}(S::SparseMatrixCSC{T,Int64},w::Vector{Float64})
  a = 0.5
  b = 1
  li = []
  lj = []
  alpha = b/(a+b)
  rtype = 2
  tol = 1e-12
  maxit = 100
  verbose = true
  ss = sum(S,2)
  P = sparse(S./ss)
  return isorank(S,w,a,b,li,lj,alpha,rtype,tol,maxit,verbose,P)
end