# Modelos de viagem espacial

Dois protótipos didáticos de pesquisa em viagem espacial, em C++ padrão (sem bibliotecas externas), com visualização no navegador:

| Módulo | Pergunta que responde | Método |
|---|---|---|
| `radiacao/` | Quanta dose de radiação uma tripulação acumula, em função da duração da missão e da espessura/material da blindagem? | Modelo semiempírico de atenuação (GCR + evento solar) |
| `reentrada/` | Como o calor de uma reentrada atravessa um escudo térmico? | Equação do calor 1D: FTCS explícito e Crank–Nicolson, validados contra solução analítica |

> **Aviso:** os parâmetros físicos (coeficientes de atenuação, condutividade, fluxo de calor) são **ilustrativos**. O repositório demonstra métodos numéricos e modelagem; não serve para projetar blindagens ou escudos reais.

## Como compilar e rodar

Precisa apenas de um compilador C++11 (`g++` ou `clang++`) e do `make`.

```bash
make run        # compila e gera os CSVs na pasta saida/
make exemplos   # copia os CSVs para docs/exemplos/ (usados pela página web)
make clean
```

Sem `make`:

```bash
g++ -O2 -std=c++11 -o radiacao radiacao/radiacao.cpp
g++ -O2 -std=c++11 -o calor1d  reentrada/calor1d.cpp
./radiacao 180 1 10     # dias, nº de eventos solares, espessura de referência (g/cm²)
./calor1d 100           # nº de intervalos; opcionalmente o dt como 2º argumento
```

## Visualização

Abra `docs/index.html` no navegador e use **Abrir meu CSV** com qualquer arquivo gerado (a 1ª coluna é o eixo X). Ativando o GitHub Pages (Settings → Pages → pasta `/docs`), os botões de exemplo também funcionam online.

## Módulo 1: radiação e blindagem

Para cada material:

```
dose_GCR_diária(x) = D0 · ( f + (1 − f) · e^(−x/λ_gcr) )      [mSv/dia]
dose_SPE(x)        = n · D_spe · e^(−x/λ_spe)                  [mSv]
dose_total         = dose_GCR_diária(x) · dias + dose_SPE(x)
```

- `x`: espessura em densidade superficial (g/cm²).
- `f`: fração da dose que a blindagem não elimina (a radiação cósmica galáctica é muito penetrante; por isso a blindagem passiva tem retorno decrescente).
- Eventos de partículas solares (SPE) são bem mais atenuados pela blindagem que a radiação galáctica.
- Materiais: alumínio, polietileno, água.

**Para tornar o modelo mais rigoroso**, calibre `D0`, `f`, `λ` com dados publicados (por exemplo as medições do instrumento RAD do Curiosity durante a viagem a Marte, Zeitlin et al., *Science*, 2013) e compare com simulações de transporte (Geant4, HZETRN).

## Módulo 2: escudo térmico de reentrada

Resolve `∂T/∂t = α ∂²T/∂x²` em `0 ≤ x ≤ L` com fluxo `q` imposto em `x = 0` e face interna isolada.

- **Explícito (FTCS):** estável apenas se `r = α·Δt/Δx² ≤ 0,5`. Rode `./calor1d 100 3` para ver a instabilidade (`r = 3`).
- **Crank–Nicolson:** estável para qualquer `Δt`, resolvido com o algoritmo de Thomas.
- **Validação:** solução analítica do meio semi-infinito com fluxo constante (Carslaw & Jaeger). Para os parâmetros padrão, o erro máximo no perfil final fica abaixo de 0,5 K. O erro na superfície é maior nos primeiros instantes porque o fluxo imposto cria uma singularidade em t → 0 (comportamento em √t).

**Próximos passos possíveis:** propriedades dependentes da temperatura, ablação, versão 2D reaproveitando a malha e a exportação VTK do solver de Poisson, estudo de convergência (erro × Δx).

## Estrutura

```
radiacao/radiacao.cpp
reentrada/calor1d.cpp
docs/index.html        # visualizador (HTML + JS puro)
docs/exemplos/         # CSVs de exemplo
docs/teoria.md         # resumo da teoria
Makefile  LICENSE  README.md
```

## Referências

- Zeitlin, C. et al. Measurements of Energetic Particle Radiation in Transit to Mars on the Mars Science Laboratory. *Science*, 340, 2013.
- Carslaw, H. S.; Jaeger, J. C. *Conduction of Heat in Solids*. Oxford University Press.
- Incropera, F. P. et al. *Fundamentos de Transferência de Calor e de Massa*. LTC.

## Licença

MIT. Veja `LICENSE`.
