abstract type KrylovSolver{T,S} end

mutable struct CgSolver{T,S} <: KrylovSolver{T,S}
  x :: S
  r :: S
  p :: S
end

function CgSolver(A, b)
  n, m = size(A)
  S = typeof(b)
  T = eltype(b)
  x = S(undef, n)
  r = S(undef, n)
  p = S(undef, n)
  return CgSolver{T,S}(x, r, p)
end

mutable struct CrSolver{T,S} <: KrylovSolver{T,S}
  x :: S
  r :: S
  p :: S
  q :: S
end

function CrSolver(A, b)
  n, m = size(A)
  S = typeof(b)
  T = eltype(b)
  x = S(undef, n)
  r = S(undef, n)
  p = S(undef, n)
  q = S(undef, n)
  return CrSolver{T,S}(x, r, p, q)
end

mutable struct SymmlqSolver{T,S} <: KrylovSolver{T,S}
  x_lq :: S
  vold :: S
  v    :: S
  w̅    :: S
end

function SymmlqSolver(A, b)
  n, m = size(A)
  S = typeof(b)
  T = eltype(b)
  x_lq = S(undef, n)
  vold = S(undef, n)
  v    = S(undef, n)
  w̅    = S(undef, n)
  return SymmlqSolver{T,S}(x_lq, vold, v, w̅)
end

mutable struct CgLanczosSolver{T,S} <: KrylovSolver{T,S}
  x      :: S
  v      :: S
  v_prev :: S
  p      :: S
end

function CgLanczosSolver(A, b)
  n, m = size(A)
  S = typeof(b)
  T = eltype(b)
  x      = S(undef, n)
  v      = S(undef, n)
  v_prev = S(undef, n)
  p      = S(undef, n)
  return CgLanczosSolver{T,S}(x, v, v_prev, p)
end

mutable struct CgLanczosShiftSolver{T,S} <: KrylovSolver{T,S}
  v          :: S
  v_prev     :: S
  x          :: Vector{S}
  p          :: Vector{S}
  σ          :: Vector{T}
  δhat       :: Vector{T}
  ω          :: Vector{T}
  γ          :: Vector{T}
  rNorms     :: Vector{T}
  indefinite :: BitArray
  converged  :: BitArray
end

function CgLanczosShiftSolver(A, b, nshifts)
  n, m = size(A)
  S = typeof(b)
  T = eltype(b)
  v          = S(undef, n)
  v_prev     = S(undef, n)
  x          = [S(undef, n) for i = 1 : nshifts]
  p          = [S(undef, n) for i = 1 : nshifts]
  σ          = Vector{T}(undef, nshifts)
  δhat       = Vector{T}(undef, nshifts)
  ω          = Vector{T}(undef, nshifts)
  γ          = Vector{T}(undef, nshifts)
  rNorms     = Vector{T}(undef, nshifts)
  indefinite = BitArray(undef, nshifts)
  converged  = BitArray(undef, nshifts)
  return CgLanczosShiftSolver{T,S}(v, v_prev, x, p, σ, δhat, ω, γ, rNorms, indefinite, converged)
end

mutable struct MinresSolver{T,S} <: KrylovSolver{T,S}
  x  :: S
  r1 :: S
  r2 :: S
  w1 :: S
  w2 :: S
end

function MinresSolver(A, b)
  n, m = size(A)
  S = typeof(b)
  T = eltype(b)
  x  = S(undef, n)
  r1 = S(undef, n)
  r2 = S(undef, n)
  w1 = S(undef, n)
  w2 = S(undef, n)
  return MinresSolver{T,S}(x, r1, r2, w1, w2)
end

mutable struct MinresQlpSolver{T,S} <: KrylovSolver{T,S}
  wₖ₋₁ :: S
  wₖ   :: S
  vₖ₋₁ :: S
  vₖ   :: S
  x    :: S
end

function MinresQlpSolver(A, b)
  n, m = size(A)
  S = typeof(b)
  T = eltype(b)
  wₖ₋₁ = S(undef, n)
  wₖ   = S(undef, n)
  vₖ₋₁ = S(undef, n)
  vₖ   = S(undef, n)
  x    = S(undef, n)
  return MinresQlpSolver{T,S}(wₖ₋₁, wₖ, vₖ₋₁, vₖ, x)
end

mutable struct DqgmresSolver{T,S} <: KrylovSolver{T,S}
  x :: S
  P :: Vector{S}
  V :: Vector{S}
  c :: Vector{T}
  s :: Vector{T}
  H :: Vector{T}
end

function DqgmresSolver(A, b, mem :: Integer=20)
  n, m = size(A)
  S = typeof(b)
  T = eltype(b)
  x = S(undef, n)
  P = [S(undef, n) for i = 1 : mem]
  V = [S(undef, n) for i = 1 : mem]
  c = Vector{T}(undef, mem)
  s = Vector{T}(undef, mem)
  H = Vector{T}(undef, mem+2)
  return DqgmresSolver{T,S}(x, P, V, c, s, H)
end

mutable struct DiomSolver{T,S} <: KrylovSolver{T,S}
  x     :: S
  x_old :: S
  P     :: Vector{S}
  V     :: Vector{S}
  L     :: Vector{T}
  H     :: Vector{T}
  p     :: BitArray
end

function DiomSolverSolver(A, b, mem :: Integer=20)
  n, m = size(A)
  S = typeof(b)
  T = eltype(b)
  x     = S(undef, n)
  x_old = S(undef, n)
  P     = [S(undef, n) for i = 1 : mem]
  V     = [S(undef, n) for i = 1 : mem]
  L     = Vector{T}(undef, mem)
  H     = Vector{T}(undef, mem+2)
  p     = BitArray(undef, mem)
  return DiomSolverSolver{T,S}(x, x_old, P, V, L, H, p)
end

mutable struct UsymlqSolver{T,S} <: KrylovSolver{T,S}
  uₖ₋₁ :: S
  uₖ   :: S
  x    :: S
  d̅    :: S
  vₖ₋₁ :: S
  vₖ   :: S
end

function UsymlqSolver(A, b)
  n, m = size(A)
  S = typeof(b)
  T = eltype(b)
  uₖ₋₁ = S(undef, m)
  uₖ   = S(undef, m)
  x    = S(undef, m)
  d̅    = S(undef, m)
  vₖ₋₁ = S(undef, n)
  vₖ   = S(undef, n)
  return UsymlqSolver{T,S}(uₖ₋₁, uₖ, x, d̅, vₖ₋₁, vₖ)
end

mutable struct UsymqrSolver{T,S} <: KrylovSolver{T,S}
  vₖ₋₁ :: S
  vₖ   :: S
  x    :: S
  wₖ₋₁ :: S
  wₖ   :: S
  uₖ₋₁ :: S
  uₖ   :: S
end

function UsymqrSolver(A, b)
  n, m = size(A)
  S = typeof(b)
  T = eltype(b)
  vₖ₋₁ = S(undef, m)
  vₖ   = S(undef, m)
  x    = S(undef, m)
  wₖ₋₁ = S(undef, m)
  wₖ   = S(undef, m)
  uₖ₋₁ = S(undef, n)
  uₖ   = S(undef, n)
  return UsymqrSolver{T,S}(vₖ₋₁, vₖ, x, wₖ₋₁, wₖ, uₖ₋₁, uₖ)
end

mutable struct TricgSolver{T,S} <: KrylovSolver{T,S}
  yₖ     :: S
  uₖ₋₁   :: S
  uₖ     :: S
  gy₂ₖ₋₁ :: S
  gy₂ₖ   :: S
  xₖ     :: S
  vₖ₋₁   :: S
  vₖ     :: S
  gx₂ₖ₋₁ :: S
  gx₂ₖ   :: S
end

function TricgSolver(A, b)
  n, m = size(A)
  S = typeof(b)
  T = eltype(b)
  yₖ     = S(undef, m)
  uₖ₋₁   = S(undef, m)
  uₖ     = S(undef, m)
  gy₂ₖ₋₁ = S(undef, m)
  gy₂ₖ   = S(undef, m)
  xₖ     = S(undef, n)
  vₖ₋₁   = S(undef, n)
  vₖ     = S(undef, n)
  gx₂ₖ₋₁ = S(undef, n)
  gx₂ₖ   = S(undef, n)
  return TricgSolver{T,S}(yₖ, uₖ₋₁, uₖ, gy₂ₖ₋₁, gy₂ₖ, xₖ, vₖ₋₁, vₖ, gx₂ₖ₋₁, gx₂ₖ)
end

mutable struct TrimrSolver{T,S} <: KrylovSolver{T,S}
  yₖ     :: S
  uₖ₋₁   :: S
  uₖ     :: S
  gy₂ₖ₋₃ :: S
  gy₂ₖ₋₂ :: S
  gy₂ₖ₋₁ :: S
  gy₂ₖ   :: S
  xₖ     :: S
  vₖ₋₁   :: S
  vₖ     :: S
  gx₂ₖ₋₃ :: S
  gx₂ₖ₋₂ :: S
  gx₂ₖ₋₁ :: S
  gx₂ₖ   :: S
end

function TrimrSolver(A, b)
  n, m = size(A)
  S = typeof(b)
  T = eltype(b)
  yₖ     = S(undef, m)
  uₖ₋₁   = S(undef, m)
  uₖ     = S(undef, m)
  gy₂ₖ₋₃ = S(undef, m)
  gy₂ₖ₋₂ = S(undef, m)
  gy₂ₖ₋₁ = S(undef, m)
  gy₂ₖ   = S(undef, m)
  xₖ     = S(undef, n)
  vₖ₋₁   = S(undef, n)
  vₖ     = S(undef, n)
  gx₂ₖ₋₃ = S(undef, n)
  gx₂ₖ₋₂ = S(undef, n)
  gx₂ₖ₋₁ = S(undef, n)
  gx₂ₖ   = S(undef, n)
  return TrimrSolver{T,S}(yₖ, uₖ₋₁, uₖ, gy₂ₖ₋₃, gy₂ₖ₋₂, gy₂ₖ₋₁, gy₂ₖ, xₖ, vₖ₋₁, vₖ, gx₂ₖ₋₃, gx₂ₖ₋₂, gx₂ₖ₋₁, gx₂ₖ)
end

mutable struct TrilqrSolver{T,S} <: KrylovSolver{T,S}
  uₖ₋₁ :: S
  uₖ   :: S
  vₖ₋₁ :: S
  vₖ   :: S
  x    :: S
  t    :: S
  d̅    :: S
  wₖ₋₁ :: S
  wₖ   :: S
end

function TrilqrSolver(A, b)
  n, m = size(A)
  S = typeof(b)
  T = eltype(b)
  uₖ₋₁ = S(undef, m)
  uₖ   = S(undef, m)
  vₖ₋₁ = S(undef, n)
  vₖ   = S(undef, n)
  x    = S(undef, m)
  t    = S(undef, n)
  d̅    = S(undef, m)
  wₖ₋₁ = S(undef, n)
  wₖ   = S(undef, n)
  return TrilqrSolver{T,S}(uₖ₋₁, uₖ, vₖ₋₁, vₖ, x, t, d̅, wₖ₋₁, wₖ)
end

mutable struct CgsSolver{T,S} <: KrylovSolver{T,S}
  x :: S
  r :: S
  u :: S
  p :: S
  q :: S
end

function CgsSolver(A, b)
  n, m = size(A)
  S = typeof(b)
  T = eltype(b)
  x = S(undef, n)
  r = S(undef, n)
  u = S(undef, n)
  p = S(undef, n)
  q = S(undef, n) 
  return CgsSolver{T,S}(x, r, u, p, q)
end

mutable struct BicgstabSolver{T,S} <: KrylovSolver{T,S}
  x :: S
  r :: S
  p :: S
  v :: S
  s :: S
end

function BicgstabSolver(A, b)
  n, m = size(A)
  S = typeof(b)
  T = eltype(b)
  x = S(undef, n)
  r = S(undef, n)
  p = S(undef, n)
  v = S(undef, n)
  s = S(undef, n)
  return BicgstabSolver{T,S}(x, r, p, v, s)
end

mutable struct BilqSolver{T,S} <: KrylovSolver{T,S}
  uₖ₋₁ :: S
  uₖ   :: S
  vₖ₋₁ :: S
  vₖ   :: S
  x    :: S
  d̅    :: S
end

function BilqSolver(A, b)
  n, m = size(A)
  S = typeof(b)
  T = eltype(b)
  uₖ₋₁ = S(undef, n)
  uₖ   = S(undef, n)
  vₖ₋₁ = S(undef, n)
  vₖ   = S(undef, n)
  x    = S(undef, n)
  d̅    = S(undef, n)
  return BilqSolver{T,S}(uₖ₋₁, uₖ, vₖ₋₁, vₖ, x, d̅)
end

mutable struct QmrSolver{T,S} <: KrylovSolver{T,S}
  uₖ₋₁ :: S
  uₖ   :: S
  vₖ₋₁ :: S
  vₖ   :: S
  x    :: S
  wₖ₋₁ :: S
  wₖ   :: S
end

function QmrSolver(A, b)
  n, m = size(A)
  S = typeof(b)
  T = eltype(b)
  uₖ₋₁ = S(undef, n)
  uₖ   = S(undef, n)
  vₖ₋₁ = S(undef, n)
  vₖ   = S(undef, n)
  x    = S(undef, n)
  wₖ₋₁ = S(undef, n)
  wₖ   = S(undef, n)
  return QmrSolver{T,S}(uₖ₋₁, uₖ, vₖ₋₁, vₖ, x, wₖ₋₁, wₖ)
end

mutable struct BilqrSolver{T,S} <: KrylovSolver{T,S}
  uₖ₋₁ :: S
  uₖ   :: S
  vₖ₋₁ :: S
  vₖ   :: S
  x    :: S
  t    :: S
  d̅    :: S
  wₖ₋₁ :: S
  wₖ   :: S
end

function BilqrSolver(A, b)
  n, m = size(A)
  S = typeof(b)
  T = eltype(b)
  uₖ₋₁ = S(undef, n)
  uₖ   = S(undef, n)
  vₖ₋₁ = S(undef, n)
  vₖ   = S(undef, n)
  x    = S(undef, n)
  t    = S(undef, n)
  d̅    = S(undef, n)
  wₖ₋₁ = S(undef, n)
  wₖ   = S(undef, n)
  return BilqrSolver{T,S}(uₖ₋₁, uₖ, vₖ₋₁, vₖ, x, t, d̅, wₖ₋₁, wₖ)
end

mutable struct CglsSolver{T,S} <: KrylovSolver{T,S}
  x :: S
  p :: S
  r :: S
end

function CglsSolver(A, b)
  n, m = size(A)
  S = typeof(b)
  T = eltype(b)
  x = S(undef, m)
  p = S(undef, m)
  r = S(undef, n)
  return CglsSolver{T,S}(x, p, r)
end

mutable struct CrlsSolver{T,S} <: KrylovSolver{T,S}
  x  :: S
  p  :: S
  Ar :: S
  r  :: S
  Ap :: S
end

function CrlsSolver(A, b)
  n, m = size(A)
  S = typeof(b)
  T = eltype(b)
  x  = S(undef, m)
  p  = S(undef, m)
  Ar = S(undef, m)
  r  = S(undef, n)
  Ap = S(undef, n)
  return CrlsSolver{T,S}(x, p, Ar, r, Ap)
end

mutable struct CgneSolver{T,S} <: KrylovSolver{T,S}
  x :: S
  p :: S
  r :: S
end

function CgneSolver(A, b)
  n, m = size(A)
  S = typeof(b)
  T = eltype(b)
  x = S(undef, m)
  p = S(undef, m)
  r = S(undef, n)
  return CgneSolver{T,S}(x, p, r)
end

mutable struct CrmrSolver{T,S} <: KrylovSolver{T,S}
  x :: S
  p :: S
  r :: S
end

function CrmrSolver(A, b)
  n, m = size(A)
  S = typeof(b)
  T = eltype(b)
  x = S(undef, m)
  p = S(undef, m)
  r = S(undef, n)
  return CrmrSolver{T,S}(x, p, r)
end

mutable struct LslqSolver{T,S} <: KrylovSolver{T,S}
  x_lq :: S
  v    :: S
  w̄    :: S
  u    :: S
end

function LslqSolver(A, b)
  n, m = size(A)
  S = typeof(b)
  T = eltype(b)
  x_lq = S(undef, m)
  v    = S(undef, m)
  w̄    = S(undef, m)
  u    = S(undef, n)
  return LslqSolver{T,S}(x_lq, v, w̄, u)
end

mutable struct LsqrSolver{T,S} <: KrylovSolver{T,S}
  x :: S
  v :: S
  w :: S
  u :: S
end

function LsqrSolver(A, b)
  n, m = size(A)
  S = typeof(b)
  T = eltype(b)
  x = S(undef, m)
  v = S(undef, m)
  w = S(undef, m)
  u = S(undef, n)
  return LsqrSolver{T,S}(x, v, w, u)
end

mutable struct LsmrSolver{T,S} <: KrylovSolver{T,S}
  x    :: S
  v    :: S
  h    :: S
  hbar :: S
  u    :: S
end

function LsmrSolver(A, b)
  n, m = size(A)
  S = typeof(b)
  T = eltype(b)
  x    = S(undef, m)
  v    = S(undef, m)
  h    = S(undef, m)
  hbar = S(undef, m)
  u    = S(undef, n)
  return LsmrSolver{T,S}(x, v, h, hbar, u)
end

mutable struct LnlqSolver{T,S} <: KrylovSolver{T,S}
  x :: S
  v :: S
  y :: S
  w̄ :: S
  u :: S
end

function LnlqSolver(A, b)
  n, m = size(A)
  S = typeof(b)
  T = eltype(b)
  x = S(undef, m)
  v = S(undef, m)
  y = S(undef, n)
  w̄ = S(undef, n)
  u = S(undef, n)
  return LnlqSolver{T,S}(x, v, y, w̄, u)
end

mutable struct CraigSolver{T,S} <: KrylovSolver{T,S}
  x :: S
  v :: S
  y :: S
  w :: S
  u :: S
end

function CraigSolver(A, b)
  n, m = size(A)
  S = typeof(b)
  T = eltype(b)
  x = S(undef, m)
  v = S(undef, m)
  y = S(undef, n)
  w = S(undef, n)
  u = S(undef, n)
  return CraigSolver{T,S}(x, v, y, w, u)
end

mutable struct CraigmrSolver{T,S} <: KrylovSolver{T,S}
  x    :: S
  v    :: S
  y    :: S
  u    :: S
  w    :: S
  wbar :: S
end

function CraigmrSolver(A, b)
  n, m = size(A)
  S = typeof(b)
  T = eltype(b)
  x    = S(undef, m)
  v    = S(undef, m)
  y    = S(undef, n)
  u    = S(undef, n)
  w    = S(undef, n)
  wbar = S(undef, n)
  return CraigmrSolver{T,S}(x, v, y, u, w, wbar)
end
