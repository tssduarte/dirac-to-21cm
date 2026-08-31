# dirac-to-21cm

Roteiros que reproduzem os resultados numéricos e a figura observacional do artigo

> Costa Neto *et al.*, *Da equação de Dirac à linha de 21 cm: um percurso didático da teoria relativística à radioastronomia*.

Repositório: <https://github.com/tssduarte/dirac-to-21cm>

Dependem apenas de `numpy`, `astropy` e `matplotlib`.

## Instalação

```bash
git clone https://github.com/tssduarte/dirac-to-21cm.git
cd dirac-to-21cm
python3 -m pip install -r requirements.txt
```

## Roteiros

| Arquivo | O que reproduz |
|---|---|
| `scripts/verifica_hfs.py` | $\Delta E_{\mathrm{hfs}}$, $A_{10}$ e a constante $\mathcal{C}$ a partir de constantes CODATA |
| `scripts/atividade_nhi.py` | atividade computacional: convenção rádio, $\int T_B\,dv$ e $N_{\mathrm{HI}}=\mathcal{C}\int T_B\,dv$ |
| `scripts/extract_hi4pi_sightlines.py` | espectros HI4PI das duas linhas de visada do artigo (área hachurada = integral; recuadro = $N_{\mathrm{HI}}$) |

```bash
python3 scripts/verifica_hfs.py
python3 scripts/atividade_nhi.py --demo
```

O modo `--demo` usa um perfil sintético e serve para auditar o encadeamento de unidades. Com um espectro real:

```bash
python3 scripts/atividade_nhi.py --fits caminho/espectro.fits
```

## Convenção de velocidade

O eixo espectral é convertido pela **convenção rádio**

$$
v = c\,(\nu_0-\nu)/\nu_0,
$$

a mesma do HI4PI, distinta da convenção óptica $v=c\,(\nu_0-\nu)/\nu$. Cubos CAR entregam `CTYPE3 = VRAD` em $\mathrm{m\,s^{-1}}$; o fator $10^3$ reduz a $\mathrm{km\,s^{-1}}$ antes de integrar. São os dois pontos em que o artigo localiza os erros mais frequentes.

## Figura HI4PI

Os cubos não entram no repositório (~250 MB cada). São produtos públicos do levantamento HI4PI ([CDS J/A+A/594/A116](https://cdsarc.cds.unistra.fr/viz-bin/qcat?J/A+A/594/A116), pasta `CUBES/EQ2000`):

- `CAR_D14.fits` — contém a linha de visada do plano galáctico
- `CAR_D09.fits` — contém a linha de visada a latitude intermediária

Coloque-os em `data/` (ou defina `HI4PI_DATA_DIR` / `--data-dir`) e rode

```bash
python3 scripts/extract_hi4pi_sightlines.py
```

Saídas: `figuras/espectros_hi4pi.pdf` e `figuras/hi4pi_numeros.tex`. Os valores publicados estão em `figuras/hi4pi_numeros.tex` e servem de referência para conferir a reprodução.

## Arquivo permanente

Uma versão etiquetada deste repositório será depositada no Zenodo, de modo a obter um DOI permanente. Até lá, cite o URL do GitHub e o artigo.
