## Desafios

### Desafio 1 - Competência básica para bioinformática

#### 1) Trimagem dos dados, por qualidade

Importar os dados

```
qiime tools import \
  --type "SampleData[SequencesWithQuality]" \
  --input-format SingleEndFastqManifestPhred33V2 \
  --input-path ./manifest.tsv \
  --output-path ./demux_seqs.qza
  ```

Visualizar os dados importados

```
qiime demux summarize \
  --i-data ./demux_seqs.qza \
  --o-visualization ./demux_seqs.qzv
  ```

A visualização mostrou gráfico (box plot) de qualidade com excelente pontuação (quality_plot.png). A parte mais baixa da box da posição 251 tinha qs=13 (apesar da mediana ser qs=37), por isso, inicialmente considerei trimar apenas a base 251. Entretanto, é bastante recomendado que se faça a trimagem de pelo menos os últimos 10 nt (que geralmente tem qualidade relativamente menor no sequenciamento Illumina). 

OBS: o relatório de qualidade completo pode ser visualizado em view.qiime2.org, usando o arquivo demux_seqs.qzv

Denoising - trim/truncate (com DADA2): DADA2 gera a feature table e as representative sequences como output, além das stats de redução de ruído (detalhadas na próxima etapa).

Sobre o dada2: 

No QIIME2, dada2 é um plugin que usa a biblioteca DADA2 do R para controle de qualidade das sequências.
O método usado foi `denoise-single` - ele toma as sequências pós demultiplexação como input; faz a redução de ruídos das sequências, realiza *dereplication* (combina as reads iguais em "sequências únicas" e fornece a quantidade de reads de cada uma dessas sequências únicas) e filtra quimeras.

- São usados os seguintes parâmetros:

--p-trunc-len: posição onde as sequências serão truncadas na extremidade 3', de acordo com a redução de qualidade (seriam as bases que foram sequenciadas nos últimos ciclos). Reads menores que o valor indicado nesse parâmetro serão descartadas.

--p-max-ee: Define o número máximo de erros esperados e descarta as reads que ultrapassem o valor. Foi usado default (2.0)

 --p-trunc-q: As reads são truncadas no quality score menor ou igual ao valor definido neste parâmetro. Se depois da truncagem a read for menor do que o definido em --p-trunc-len, ela é descartada. Foi usado default (2).

 --p-pooling-method: Determina o método usado para fazer o pool das amostras para redução de ruído. Pode ser escolhida entre as opções "independent" e "pseudo". Foi usado independent (default), que faz redução de ruído das amostras de forma independente.

 --p-chimera-method: Determina o método para remoção de quimeras, podendo ser "pooled", "none" ou "consensus". Foi usado o default, consensus, no qual as quimeras são detectadas nas amostras individualmente e são removidas as sequências quiméricas que estejam presentes em uma quantidade suficiente de amostras.

 --p-min-fold-parent-over-abundance: Este parâmetro determina a abundância mínima que devem ter as sequências parentais de uma sequência que está sendo testada como quimérica. A medida é expressada como a razão entre a sequência parental e a sequência que está sendo testada e seu valor deve ser maior ou igual a 1, ou seja, as sequências parentais devem ser mais abundantes que as sequências testadas como quiméricas. Foi usado default (1).

 --p-n-reads-learn: Determina o número de reads a ser usado para treinar o modelo de erro. Foi usado default (1000000).

- São gerados como outputs:

--o-table: artefato do qiime relativo à feature table resultante

--o-representative-sequences: artefato do qiime relativo às sequências que representam as features. Cada feature da feature table será representada por uma única sequência.

--o-denoising-stats: relatório de redução de ruído, detalhado abaixo.

 
  ```
  qiime dada2 denoise-single \
  --i-demultiplexed-seqs ./demux_seqs.qza \
  --p-trunc-len 241 \                        
  --o-table ./dada2_table241.qza \
  --o-representative-sequences ./dada2_rep_set241.qza \
  --o-denoising-stats ./dada2_stats241.qza
  ```

Visualizar o relatório e estatísticas do denoising:

Usando o arquivo dada2_stats241.qzv em view.qiime2.org é possível ver, para cada amostra, a quantidade de sequências inseridas como input, a quantidade restante após filtragem (e sua porcentagem), a quantidade de sequências que passaram por redução de ruído e a quantidade e porcentagem de input que não continha quimeras.

```
qiime metadata tabulate \
  --m-input-file ./dada2_stats241.qza  \
  --o-visualization ./dada2_stats241.qzv
```

Gerar uma versão tabulada da feature table para visualização:

```
qiime metadata tabulate \
  --m-input-file dada2_table241.qza \
  --o-visualization feature_table.qzv
```

#### 3) Identificação taxonômica

Usando classificador treinado Greengenes 13_8 99% OTUs da região 515F/806R (genes de 16s clusterizados a 99% de semelhança, com o seguinte par de primers: forward a partir da posição 515 e reverse na posição 806, ou seja, abrangem a região V4 do gene 16s). Nessa etapa, QIIME2 trabalha atrvés da scikit-learn, que pede uso do classificador treinado, já que se trata de uma biblioteca de machine learning.

OBS: Tentei usar o SILVA, mas a execução do código não se completou, mostrando apenas a mensagem "killed" - acredito estar ligada à limitação da capacidade de processamento da máquina. Por isso usei o Greengenes treinado, mais leve.

  ```
  qiime feature-classifier classify-sklearn \
  --i-reads ./dada2_rep_set241.qza \
  --i-classifier ./gg-13-8-99-515-806-nb-classifier.qza \
  --o-classification ./taxonomy.qza
  ```
  
  Visualização da identificação taxonômica
  Em view.qiime2.org, usando o arquivo taxonomy.qzv

  ```
  qiime metadata tabulate \
  --m-input-file ./taxonomy.qza \
  --o-visualization ./taxonomy.qzv
  ```

O método a seguir gera um gráfico de barras interativo que mostra a identificação taxonômica em cada amostra (visualizar o arquivo QZV em view.qiime2.org).

  ```
  qiime taxa barplot \
  --i-table dada2_table241.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file cmetadata.tsv \
  --o-visualization taxa-bar-plots.qzv
  ```


#### 4) Gerar a OTU table com linhas sendo taxonomias e colunas sendo as amostras

Colapsar as features: combinar as sequências que tenham a mesma taxonomia, ainda que sejam ligeiramente diferentes. Nesse caso, as frequências das features são somadas.

```
qiime taxa collapse \
    --i-table table.qza \           #input: feature table gerada na etapa do dada2 (com o nome alterado)
    --i-taxonomy taxonomy.qza \     #input: taxonomias atribuídas às features da feature table fornecida
    --p-level 7 \
    --o-collapsed-table table-l7.qza #output: feature table resultante, onde todas as features agora são identificações taxonômicas com o nível taxonômico indicado no parâmetro acima (p-level 7)
```

Transpor a feature table gerada no collapse (transformar linhas em colunas e colunas em linhas)

``` 
qiime feature-table transpose \
--i-table tablel7.qza \
--o-transposed-feature-table trans_table.qza
```

Gerar a visualização da tabela final (**final_table.tsv**)

```
qiime metadata tabulate \
    --m-input-file trans_table.qza \
    --o-visualization final_table.qzv #este output serve para visualização em view.qiime2.org da mesma tabela TSV
```







  

