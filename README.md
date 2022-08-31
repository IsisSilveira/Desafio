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

  A visualização mostrou gráfico (box plot) de qualidade com excelente pontuação. A parte mais baixa da box da posição 251 tinha qs=13 (apesar da mediana ser qs=37), por isso, inicialmente considerei trimar apenas a base 251. Entretanto, é bastante recomendado que se faça a trimagem de pelo menos os últimos 10 nt (que geralmente tem qualidade relativamente menor no sequenciamento Illumina). 
  OBS: o relatório de qualidade completo pode ser visualizado em view.qiime2.org, usando o arquivo demux_seqs.qzv

  Denoising - trim/truncate (com DADA2)
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

  Classificação taxonômica

  ```
  qiime feature-classifier classify-sklearn \
  --i-reads ./dada2_rep_set241.qza \
  --i-classifier ./gg-13-8-99-515-806-nb-classifier.qza \
  --o-classification ./taxonomy.qza
  ```
  
  Visualização da classificação taxonômica
  Em view.qiime2.org, usando o arquivo taxonomy.qzv

  ```
  qiime metadata tabulate \
  --m-input-file ./taxonomy.qza \
  --o-visualization ./taxonomy.qzv
  ```






  

