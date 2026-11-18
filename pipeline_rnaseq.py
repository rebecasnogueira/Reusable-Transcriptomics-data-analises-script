"""Pipeline reutilizável de exploração de matrizes de contagens RNA-seq.

Instalação: python -m pip install pandas numpy matplotlib seaborn scikit-learn
Uso: python pipeline_rnaseq_reutilizavel.py --rna contagens.txt
     python pipeline_rnaseq_reutilizavel.py --rna contagens.csv --metadata metadata.csv
Ajuda: python pipeline_rnaseq_reutilizavel.py --help

Entrada: genes nas linhas, amostras nas colunas, contagens brutas não negativas.
Formatos: CSV/TSV com cabeçalho padrão, ou TXT separado por espaços com
cabeçalho sem nome para a primeira coluna (formato do rawcountAbate.txt).
Metadados opcionais: colunas de ID da amostra e grupo, configuráveis.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA

# ============ CONFIGURAÇÕES PADRÃO (EDITE AQUI OU USE --help) ============
ARQUIVO_RNA = 'rawcountAbate.txt'
ARQUIVO_METADATA = None  # Exemplo: 'metadata.csv'
PASTA_SAIDA = 'Resultados_RNAseq'
SEPARADOR_RNA = 'auto'   # auto | virgula | tab | espacos
CABECALHO_SEM_ID = False  # True para rawcountAbate.txt (ou auto detectado)
CONTAGEM_MINIMA = 10
MINIMO_AMOSTRAS = 3
COLUNA_AMOSTRA = '#NAME'
COLUNA_GRUPO = 'class'
METODO_CORRELACAO = 'pearson'
DPI = 250
# =======================================================================


def importar_contagens(caminho: Path, separador: str = 'auto', cabecalho_sem_id: bool = False) -> pd.DataFrame:
    if not caminho.is_file():
        raise FileNotFoundError(f'Arquivo não encontrado: {caminho}')
    sep = {'auto': None, 'virgula': ',', 'tab': '\t', 'espacos': r'\s+'}[separador]
    if sep is None:
        if caminho.suffix.lower() == '.csv':
            sep = ','
        elif caminho.suffix.lower() in {'.tsv', '.tab'}:
            sep = '\t'
        else:
            sep = r'\s+'
    if sep == r'\s+' and (cabecalho_sem_id or separador == 'auto'):
        with caminho.open(encoding='utf-8-sig') as f:
            header = f.readline().split()
            primeira_linha = f.readline().split()
        # Detecta cabeçalho sem coluna de gene: dados possuem uma coluna a mais.
        if len(primeira_linha) == len(header) + 1:
            rna = pd.read_csv(caminho, sep=r'\s+', skiprows=1, header=None, index_col=0)
            rna.columns = header
        else:
            rna = pd.read_csv(caminho, sep=sep, index_col=0)
    else:
        rna = pd.read_csv(caminho, sep=sep, index_col=0)
    rna.index = rna.index.astype(str).str.strip()
    rna.columns = rna.columns.astype(str).str.strip()
    rna.index.name = 'Gene_ID'
    return rna


def controle_qualidade(rna: pd.DataFrame, saida: Path) -> pd.DataFrame:
    if rna.empty or rna.shape[1] < 2:
        raise ValueError('A matriz deve conter genes e pelo menos duas amostras.')
    dados = rna.apply(pd.to_numeric, errors='coerce')
    erros = {
        'valores_ausentes_ou_nao_numericos': int(dados.isna().sum().sum()),
        'genes_duplicados': int(rna.index.duplicated().sum()),
        'amostras_duplicadas': int(rna.columns.duplicated().sum()),
        'valores_negativos': int((dados < 0).sum().sum()),
        'valores_nao_inteiros': int(((dados % 1 != 0) & dados.notna()).sum().sum()),
    }
    pd.Series(erros, name='Quantidade').to_csv(saida / 'QC_validacao.csv', header=True)
    if any(erros.values()):
        raise ValueError(f'Matriz inválida: {erros}. Verifique QC_validacao.csv')
    resumo = pd.DataFrame({
        'Total_contagens': dados.sum(),
        'Genes_detectados': (dados > 0).sum(),
        'Genes_zero': (dados == 0).sum(),
        'Media_contagens': dados.mean(),
        'Mediana_contagens': dados.median(),
    })
    resumo['Percentual_detectados'] = 100 * resumo['Genes_detectados'] / len(dados)
    resumo.index.name = 'Amostra'
    resumo.to_csv(saida / 'QC_resumo_amostras.csv')
    return resumo


def salvar(fig, saida: Path, nome: str) -> None:
    fig.savefig(saida / nome, dpi=DPI, bbox_inches='tight')
    plt.close(fig)


def grafico_bibliotecas(resumo: pd.DataFrame, saida: Path) -> None:
    fig, ax = plt.subplots(figsize=(max(9, len(resumo) * .65), 5))
    ax.bar(resumo.index, resumo['Total_contagens'] / 1e6)
    ax.set(ylabel='Contagens (milhões)', title='Tamanho das bibliotecas')
    ax.tick_params(axis='x', labelrotation=55)
    fig.tight_layout()
    salvar(fig, saida, 'QC_tamanho_bibliotecas.png')


def filtrar(rna: pd.DataFrame, contagem: int, amostras: int) -> pd.DataFrame:
    if contagem < 0 or not 1 <= amostras <= rna.shape[1]:
        raise ValueError('Parâmetros de filtragem inválidos.')
    filtrado = rna.loc[(rna >= contagem).sum(axis=1) >= amostras].copy()
    if filtrado.empty:
        raise ValueError('Nenhum gene permaneceu após a filtragem.')
    return filtrado


def normalizar(rna: pd.DataFrame, filtrado: pd.DataFrame) -> pd.DataFrame:
    totais = rna.sum(axis=0)  # Total original, antes da filtragem.
    if (totais <= 0).any():
        raise ValueError('Biblioteca com zero contagens totais.')
    return filtrado.div(totais, axis=1) * 1_000_000


def grafico_boxplot(log2cpm: pd.DataFrame, saida: Path) -> None:
    fig, ax = plt.subplots(figsize=(max(10, log2cpm.shape[1] * .7), 6))
    ax.boxplot([log2cpm[c].to_numpy() for c in log2cpm], tick_labels=list(log2cpm), showfliers=False)
    ax.set(ylabel='log2(CPM + 1)', title='Distribuição da expressão gênica por amostra')
    ax.tick_params(axis='x', labelrotation=55)
    fig.tight_layout()
    salvar(fig, saida, 'boxplot_log2_CPM.png')


def correlacoes(log2cpm: pd.DataFrame, saida: Path, metodo: str) -> None:
    cor = log2cpm.corr(method=metodo)
    cor.to_csv(saida / 'matriz_correlacao.csv', index_label='Amostra')
    sem_diag = cor.mask(np.eye(len(cor), dtype=bool))  # Sem alteração de array read-only.
    sem_diag.mean(axis=1).rename('Correlacao_media').to_csv(saida / 'correlacao_media.csv')
    valores = sem_diag.to_numpy(copy=True)
    minimo = float(np.nanmin(valores))
    for nome, vmin in [('heatmap_correlacao_0a1.png', 0), ('heatmap_correlacao_ajustado.png', minimo)]:
        fig, ax = plt.subplots(figsize=(max(9, len(cor) * .7), max(7, len(cor) * .65)))
        sns.heatmap(cor, annot=len(cor) <= 25, fmt='.3f', cmap='viridis', vmin=vmin, vmax=1,
                    square=True, linewidths=.3, cbar_kws={'label': f'Correlação ({metodo})'}, ax=ax)
        ax.set_title('Correlação entre amostras' + (' — escala ajustada' if vmin == minimo else ''))
        ax.tick_params(axis='x', labelrotation=45)
        fig.tight_layout()
        salvar(fig, saida, nome)


def carregar_grupos(caminho: Path | None, ids: pd.Index, col_id: str, col_grupo: str) -> pd.Series | None:
    if caminho is None:
        return None
    meta = pd.read_csv(caminho)
    if col_id not in meta or col_grupo not in meta:
        raise ValueError(f'Metadados precisam conter as colunas {col_id!r} e {col_grupo!r}.')
    meta[col_id] = meta[col_id].astype(str).str.strip()
    if meta[col_id].duplicated().any():
        raise ValueError('Metadados possuem identificadores de amostras duplicados.')
    grupos = meta.set_index(col_id)[col_grupo].reindex(ids)
    if grupos.isna().any():
        print('AVISO: amostras sem grupo serão mostradas como "Sem grupo":', grupos[grupos.isna()].index.tolist())
        grupos = grupos.fillna('Sem grupo')
    return grupos.astype(str)


def executar_pca(log2cpm: pd.DataFrame, saida: Path, grupos: pd.Series | None) -> None:
    dados = log2cpm.T
    dados = dados.loc[:, dados.var(axis=0) > 0]
    if dados.shape[1] < 2 or dados.shape[0] < 3:
        print('PCA ignorada: amostras ou genes variáveis insuficientes.')
        return
    if not np.isfinite(dados.to_numpy()).all():
        raise ValueError('Dados não finitos impedem PCA.')
    modelo = PCA(n_components=2)
    coords = modelo.fit_transform(dados)
    variancia = modelo.explained_variance_ratio_ * 100
    pcs = pd.DataFrame(coords, index=dados.index, columns=['PC1', 'PC2'])
    pcs.index.name = 'Amostra'
    pcs.to_csv(saida / 'PCA_coordenadas.csv')
    pd.DataFrame({'Componente': ['PC1', 'PC2'], 'Variancia_explicada_percentual': variancia}).to_csv(
        saida / 'PCA_variancia_explicada.csv', index=False)
    if grupos is not None:
        pcs['Grupo'] = grupos.reindex(pcs.index)
        pcs.to_csv(saida / 'PCA_coordenadas_com_grupos.csv')
    fig, ax = plt.subplots(figsize=(10, 7))
    if grupos is None:
        ax.scatter(pcs.PC1, pcs.PC2, s=85)
    else:
        for grupo, bloco in pcs.groupby('Grupo'):
            ax.scatter(bloco.PC1, bloco.PC2, s=85, label=str(grupo))
        ax.legend(title='Grupo')
    for amostra, linha in pcs.iterrows():
        ax.annotate(str(amostra), (linha.PC1, linha.PC2), xytext=(4, 4), textcoords='offset points', fontsize=8)
    ax.axhline(0, color='gray', ls='--', lw=.7)
    ax.axvline(0, color='gray', ls='--', lw=.7)
    ax.set(xlabel=f'PC1 ({variancia[0]:.2f}%)', ylabel=f'PC2 ({variancia[1]:.2f}%)', title='PCA de RNA-seq')
    fig.tight_layout()
    salvar(fig, saida, 'PCA_RNAseq.png' if grupos is None else 'PCA_RNAseq_por_grupo.png')


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--rna', default=ARQUIVO_RNA, help='Matriz de contagens')
    parser.add_argument('--metadata', default=ARQUIVO_METADATA, help='CSV opcional com grupos')
    parser.add_argument('--saida', default=PASTA_SAIDA, help='Pasta de resultados')
    parser.add_argument('--separador', choices=['auto', 'virgula', 'tab', 'espacos'], default=SEPARADOR_RNA)
    parser.add_argument('--cabecalho-sem-id', action='store_true', default=CABECALHO_SEM_ID)
    parser.add_argument('--contagem-minima', type=int, default=CONTAGEM_MINIMA)
    parser.add_argument('--minimo-amostras', type=int, default=MINIMO_AMOSTRAS)
    parser.add_argument('--coluna-amostra', default=COLUNA_AMOSTRA)
    parser.add_argument('--coluna-grupo', default=COLUNA_GRUPO)
    parser.add_argument('--correlacao', choices=['pearson', 'spearman'], default=METODO_CORRELACAO)
    args = parser.parse_args()
    saida = Path(args.saida).expanduser().resolve()
    saida.mkdir(parents=True, exist_ok=True)
    rna = importar_contagens(Path(args.rna).expanduser(), args.separador, args.cabecalho_sem_id)
    resumo = controle_qualidade(rna, saida)
    grafico_bibliotecas(resumo, saida)
    filtrado = filtrar(rna, args.contagem_minima, args.minimo_amostras)
    cpm = normalizar(rna, filtrado)
    log2cpm = np.log2(cpm + 1)
    log2cpm.to_csv(saida / 'RNAseq_log2_CPM.csv', index_label='Gene_ID')
    pd.Series({'Genes_originais': len(rna), 'Genes_filtrados': len(filtrado),
               'Genes_removidos': len(rna) - len(filtrado), 'Amostras': rna.shape[1],
               'Contagem_minima': args.contagem_minima, 'Minimo_amostras': args.minimo_amostras},
              name='Valor').to_csv(saida / 'resumo_processamento.csv', header=True)
    grafico_boxplot(log2cpm, saida)
    correlacoes(log2cpm, saida, args.correlacao)
    grupos = carregar_grupos(Path(args.metadata).expanduser() if args.metadata else None,
                             log2cpm.columns, args.coluna_amostra, args.coluna_grupo)
    executar_pca(log2cpm, saida, grupos)
    print(f'Concluído: {rna.shape[0]} genes iniciais; {len(filtrado)} após filtro; {rna.shape[1]} amostras.')
    print(f'Resultados: {saida}')


if __name__ == '__main__':
    main()
