# Resumo da teoria

## Radiação no espaço profundo
Duas fontes dominam: raios cósmicos galácticos (GCR), contínuos e muito energéticos, e eventos de partículas solares (SPE), esporádicos e menos penetrantes. A blindagem por materiais leves e ricos em hidrogênio (polietileno, água) tende a funcionar melhor que metais pesados, porque gera menos partículas secundárias por fragmentação nuclear. Por isso o modelo usa uma fração residual `f` maior para o alumínio.

## Condução de calor
A equação do calor `∂T/∂t = α ∂²T/∂x²` descreve a difusão térmica. Discretizando com diferenças finitas (`Δx`, `Δt`) e `r = αΔt/Δx²`:

- **FTCS:** `T_i^{n+1} = T_i^n + r (T_{i+1}^n − 2T_i^n + T_{i−1}^n)`. Estabilidade de von Neumann exige `r ≤ 1/2`.
- **Crank–Nicolson:** média entre os níveis `n` e `n+1`; precisão de 2ª ordem no tempo e no espaço e estabilidade incondicional. Gera um sistema tridiagonal em cada passo.
- **Condição de fluxo em x = 0:** nó fantasma `T_{−1} = T_1 + 2Δx·q/k`, que dá a equação do nó de contorno com precisão de 2ª ordem.
