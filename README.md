# dirac-to-21cm

Da equação de Dirac à densidade em coluna de HI: derivação da linha de 21 cm e
aplicação a espectros públicos do HI4PI.

Roteiros que reproduzem os resultados numéricos e a figura observacional de

> Costa Neto *et al.*, *Da equação de Dirac à linha de 21 cm: um percurso
> didático da teoria relativística à radioastronomia*. Submetido à Revista
> Brasileira de Ensino de Física.

## Como citar

Cite o artigo. Para referir-se ao código em si, use também este repositório:

```bibtex
@software{dirac21cm,
  author  = {Costa Neto, ... e Duarte, T. S. S.},
  title   = {dirac-to-21cm: da equação de Dirac à densidade em coluna de HI},
  year    = {2026},
  url     = {https://github.com/tssduarte/dirac-to-21cm}
}
```

Uma versão etiquetada será depositada no Zenodo para obter DOI permanente; a
entrada acima será atualizada quando isso ocorrer.

## Roteiros

| Arquivo | O que reproduz |
|---|---|
| `scripts/verifica_hfs.py` | $\Delta E_{\mathrm{hfs}}$, $A_{10}$ e a constante $\mathcal{C}$, a partir de constantes CODATA (Tabela 1 do artigo) |
| `scripts/atividade_nhi.py` | atividade computacional: convenção rádio, $\int T_B\,dv$ e $N_{\mathrm{HI}}=\mathcal{C}\int T_B\,dv$ (Seção 7) |
| `scripts/extract_hi4pi_sightlines.py` | espectros das duas linhas de visada, Figura 5 e Tabela 2 |

## Uso

```bash
git clone https://github.com/tssduarte/dirac-to-21cm.git
cd dirac-to-21cm
python3 -m pip install -r requirements.txt

python3 scripts/verifica_hfs.py
python3 scripts/atividade_nhi.py --demo
```

`verifica_hfs.py` não depende de dados externos e reproduz a Tabela 1 em
segundos — é a verificação mais rápida de que a instalação funciona.

O modo `--demo` usa um perfil sintético e serve para auditar o encadeamento de
unidades. Com um espectro real:

```bash
python3 scripts/atividade_nhi.py --fits caminho/espectro.fits
```

## Dados do HI4PI

Os cubos não entram no repositório (~250 MB cada). São produtos públicos do
levantamento HI4PI ([CDS J/A+A/594/A116](https://cdsarc.cds.unistra.fr/viz-bin/qcat?J/A+A/594/A116),
pasta `CUBES/EQ2000`):

- `CAR_D14.fits` — linha de visada do plano galáctico
- `CAR_D09.fits` — linha de visada a latitude intermediária

Coloque-os em `data/` (ou defina `HI4PI_DATA_DIR` / `--data-dir`) e rode:

```bash
python3 scripts/extract_hi4pi_sightlines.py
```

Saídas: `figuras/espectros_hi4pi.pdf` e `figuras/hi4pi_numeros.tex`.
O PDF versionado em `figuras/espectros_hi4pi.pdf` é o gabarito da Figura 5
(área hachurada = $\int T_B\,dv$; recuadro = $N_{\mathrm{HI}}$); as macros em
`figuras/hi4pi_numeros.tex` são o gabarito numérico da Tabela 2.

## Convenção de velocidade

O eixo espectral é convertido pela **convenção rádio**,

$$v = c\,(\nu_0-\nu)/\nu_0,$$

a mesma adotada pelo HI4PI e distinta da convenção óptica
$v = c\,(\nu_0-\nu)/\nu$. Os cubos CAR entregam `CTYPE3 = VRAD` em
m s⁻¹; os roteiros convertem para km s⁻¹ antes de integrar. Uma ordem de
grandeza de $N_{\mathrm{HI}}$ fora da faixa $10^{19}$–$10^{22}$ cm⁻² costuma
indicar erro nessa conversão.

## Licença

Código sob licença MIT (ver `LICENSE`). Figuras e material didático sob
CC BY 4.0.
