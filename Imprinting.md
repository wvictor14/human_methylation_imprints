Imprinting is a biological feature where a only one allele is expressed
for a gene, which is dependent on its parental origin. This
allele-specific expression is regulated by allele-specific DNA
methylation.

Most imprinted genes are ubiquitous, the allele-specificity of
expression and DNA methylation is observed in all tissues. However, the
placenta also expresses many genes in an allele-specific parental-origin
dependent manner. This indicates that imprinting in the placenta serves
important roles to it’s function, which is not surprising, given that
with respect to evolution, imprinting and placentation evolved at the
same time.

Many studies have characterized imprinting using expression technology
(e.g. RNAseq) and DNA methylation. These studies have accumulated
evidence for the existence of many imprints. However, there is lacking a
universal reference for human imprinted genes and regions.

Here, I try to compile the evidence for all detected human imprinted
regions. My goal is to create a hg19 reference that contains the
coordinates for imprinted regions that are 1) ubiquitously expressed in
all / most tissues, and 2) specific to the placenta.

For each publication, I had to do some extra work to obtain the
necessary information on the imprinted regions. I describe this extra
processing in each section. But briefly, I had to do things like: figure
out which regions were confirmed to be allelically-methylated, figure
out which are placental-specific or not.

Primary Publications
====================

There are many studies that have looked at specific imprinted regions
and genes. Here, I briefly describe the most comprehensive studies that
have been done to date. Most of these have used some sort of
next-generation sequencing technology in combination with high density
microarrays to measure DNA methylation and sometimes transcription.

Hanna et al 2016
----------------

**Genome Build:** hg19 **Tissues used:** Placenta, Germ cells
**Technology:** 450k Array **Approach to determine parental origin:**
DMR analysis of triploid (diandric vs digynic) conceptions overlapping
with DMRs between oocyte and sperm **Findings:** Imprinted regions
identified for known genes. Novel imprinted regions too. Polymorphic
placental-specific imprinting. Some Court et al 2014 DMRs were not
reproduced, and did not show germline differences in methylation.
**Number of DMRs:** 43 somatic imprinted DMRs in the placenta that
overlapped to Court et al. 2014. 10 of which may be secondary DMRs (no
germline differences). 101 novel maternal germline imprinted DMRs, 72
are placenta specific.

At this point &gt; 140 maternal imprinted DMRs in placenta.

limitation is small sample numbers, discovery in 450k CpGs, so there are
likely many more DMRs.

> Germline novel DMRs were those that were CGIs (defined by UCSC Genome
> Browser
> \[<a href="https://genome.ucsc.edu/" class="uri">https://genome.ucsc.edu/</a>\])
> that (1) were &gt;50% differentially methylated between sperm and
> oocytes, (2) were intermediately methylated (15%–60%, based on the
> 90th percentile observed at known imprinted DMRs) in the blastocyst
> (average of ICM and TE), (3) showed a &gt;5% difference in methylation
> between diandric and digynic triploids, and (4) showed matching
> parental origin of DNA methylation between triploid villi and gametes.
> Placental-specific DMRs were defined as those with &gt;25% and &lt;75%
> DNA methylation in placenta and &lt;25% in somatic tissues.

    hanna_all <- read_excel(here::here('processed', 'Hanna_processed.xlsx'), skip = 1, sheet = 1) %>%
      janitor::clean_names()
    hanna_all

    ## # A tibble: 144 x 9
    ##      dmr chr_start_end closest_tss region hgnc_id cgi_density number_probes
    ##    <dbl> <chr>         <chr>       <chr>  <chr>   <chr>               <dbl>
    ##  1     1 chr1:4002497… PPIEL       Promo… 33195   HC                      9
    ##  2     2 chr1:6851183… DIRAS3 (DM… Promo… 687     HC                     10
    ##  3     3 chr1:6851578… DIRAS3 (DM… Promo… 687     HC                     13
    ##  4     4 chr4:8961792… NAP1L5      Promo… 19968   HC                      8
    ##  5     5 chr6:3848898… FAM50B      Promo… 18789   HC                     28
    ##  6     6 chr6:1443284… PLAGL1      Promo… 9046    HC                     18
    ##  7     7 chr7:9428550… PEG10/SGCE  Promo… 14005/… HC                     55
    ##  8     8 chr7:1301303… MEST/MESTI… Promo… 7028/1… HC                     56
    ##  9     9 chr7:1548615… HTR5A       Promo… 5300    HC                     14
    ## 10    10 chr8:3760551… ERLIN2      Promo… 1356    IC                      6
    ## # … with 134 more rows, and 2 more variables: methylated_allele <chr>,
    ## #   court_et_al_2014_dmr <chr>

    hanna_pl <- read_excel(here::here('processed', 'Hanna_processed.xlsx'), skip = 1, sheet = 2) %>%
      janitor::clean_names()
    hanna_pl

    ## # A tibble: 88 x 13
    ##    gene_name ucsc_coordinate… chromosome  start    end dmr_id ensembl_gene
    ##    <chr>     <chr>            <chr>       <dbl>  <dbl> <chr>  <chr>       
    ##  1 THAP3     chr1:6685210-66… chr1       6.69e6 6.69e6 THAP3  THAP3       
    ##  2 AKR7A3    chr1:19609057-1… chr1       1.96e7 1.96e7 AKR7A3 AKR7A3      
    ##  3 C1orf216  chr1:36179477-3… chr1       3.62e7 3.62e7 C1orf… C1orf216    
    ##  4 ACOT11    chr1:55013807-5… chr1       5.50e7 5.51e7 ACOT11 ACOT11      
    ##  5 PCSK9     chr1:55505149-5… chr1       5.55e7 5.55e7 PCSK9  PCSK9       
    ##  6 IL12RB2   chr1:67773047-6… chr1       6.78e7 6.79e7 IL12R… IL12RB2     
    ##  7 TNR       chr1:175291935-… chr1       1.75e8 1.76e8 TNR    TNR         
    ##  8 CACNA1E   chr1:181452686-… chr1       1.81e8 1.82e8 CACNA… CACNA1E     
    ##  9 G0S2      chr1:209848670-… chr1       2.10e8 2.10e8 G0S2   G0S2        
    ## 10 LINC00467 chr1:211556097-… chr1       2.12e8 2.12e8 LINC0… LINC00467   
    ## # … with 78 more rows, and 6 more variables: ensembl_chr <dbl>,
    ## #   ensembl_start <dbl>, ensembl_end <dbl>, te <chr>,
    ## #   x2nd_trimester_placenta <chr>, x3rd_trimester_placenta <chr>

    hanna_all %>%
      mutate(placental_specific = closest_tss %in% hanna_pl$dmr_id,
             known_dmr = court_et_al_2014_dmr != 'NA') %>%
      count(known_dmr, placental_specific)

    ## # A tibble: 4 x 3
    ##   known_dmr placental_specific     n
    ##   <lgl>     <lgl>              <int>
    ## 1 FALSE     FALSE                 29
    ## 2 FALSE     TRUE                  72
    ## 3 TRUE      FALSE                 27
    ## 4 TRUE      TRUE                  16

    hanna_final <- hanna_all %>%
      mutate(placental_specific = closest_tss %in% hanna_pl$dmr_id,
             known_dmr = court_et_al_2014_dmr != 'NA') 

    hanna_final[hanna_final == 'NA'] <- NA

Note, of 101 novel imprinted DMRs, 72 are placental specific. These are
not indicated in s2, but were used for the DAVID analysis in s6. Also,
coordinates from s6 do not correspond exactly to s2. In most cases,
there are a couple 100 bp differences maximum. I’m not sure exactly but
I suspect that s6 is the coordinates of the main isoform/transcript
instead of the DMR. I use all coordinates from s2.

    hanna_final %>%
      write_tsv(here::here('processed', 'hanna_final.tsv'))

**definition of placental-specific:** &gt; Seventy-two novel DMRs showed
a placental-specific pattern of imprinting, defined as &lt;25% or
&gt;75% methylation in somatic tissues and intermediately (25%–75%)
methylated in placenta.

Where somatic tissues are fetal brain, spinal cord, kidney, muscle, and
adult blood.

Court et al 2014
----------------

**Tissues:** Placenta, Normal human tissues **Technology:** WGBSEQ, 450k
**Approach to determine parental origin:** UPDs and hydatidiform moles.
**Findings:** Imprinted regions for all known imprinted genes. 21 novel
loci.

Limitation discovery in 450k CpGs (intersected with wgBSeq).

It appears that all relevant data on the confirmed regions is in table
1, and not in the supplementals, which are terribly terribly annotated.

    court_ubi <- read_excel(here::here('processed','Court_processed.xlsx'), sheet = 1)
    court_pl <- read_excel(here::here('processed','Court_processed.xlsx'), sheet = 2)

    court <- court_ubi %>% mutate(type = 'ubiquitous') %>%
      bind_rows(court_pl %>% mutate(type = 'placental-specific')) %>%
      janitor::clean_names() %>%
      select(chr, start, finish, type, known_dmr, known_imprinting_loci, everything()) %>%
      dplyr::rename(end = finish, number_450k_probes = number_infinium_probes, 
                    methylated_allele = methylation_origin)

    court %>%
      write_tsv(here::here('processed', 'court_final.tsv'))

Sanchez-Delgado et al 2016
--------------------------

**Genome build:** hg19 **Tissues:** germ cells, blastocysts, placenta,
various somatic tissues **Technology:** Methylseq **Approach:**
Identified intermediately methylated PL regions, that are also DMRs in
germ cells. Validated with allele-specific assays in 11% of 551
identified regions. **Findings:** Sometimes imprinted DMRs do not result
in imprinted expression, confirm derived imprinted methylation is
maternal/pat with usage of germ cells methylation data. **Number of
DMRs:** Verified Hanna and Court 150 maternally methylated
placental-specific DMRs (previously idetnfieid). 551 additional novel
DMRs candidates. Only 11% were verified to be monoallelically methylated
(all maternal).

Table S2 has a list of regions. They list all candidate 551 maternal
gDMRs. 11% of these were confirmed to be mono-allelic methylated. These
are not indicated in table S2, but are indicated in table S3. I need
table S2 for coordinates of the DMRs.

Table S3 needs to get processed further. I notice that the sample
numbers are really small. It’s not convincing that these DMRs are
imprinted.

We filter to DMRs with “enough” evidence to support their status:

-   number of informative cases is greater or equal to 3

After that, we can define impritning as,

-   maternall methylated when \# of maternal = n\_all
-   polymorphic when n maternally methylated &lt; n\_all

**Note** that they do not provide coordinates for the verified DMRs in
table S3. Therefore, I have to merge S2 and S3 based on the only column
in common, “Gene”. And note that the “Gene” column in table S3 contains
non-specific information such as “Chr18” and “Chr1”???. As expected,
many of these cannot be merged into table S2 and will just be dropped
from the analysis..

    sandel_candidates <- read_excel(here::here('processed', 'Sanchez_Delgado_processed.xlsx'), sheet = 1)
    sandel_confirmed <- read_excel(here::here('processed', 'Sanchez_Delgado_processed.xlsx'), sheet = 2)

    sandel_confirmed <- sandel_confirmed %>% 
      # extract evidence for imprinting
      dplyr::rename(description = names(.)[2]) %>%
      mutate(n_mat_methylated = str_extract(description, '[0-9](?= mat)') %>% as.numeric(),
             n_pat_methylated =  str_extract(description, '[0-9](?= pat)') %>% as.numeric(),
             n_bi_methylated =  str_extract(description, '[0-9](?= biall)') %>% as.numeric()) %>%
      rowwise() %>%
      mutate(n_all = sum(n_mat_methylated, n_pat_methylated, n_bi_methylated, na.rm = T)) %>%
      ungroup() %>%
      
      # filter to where there are enough cases to make an inference about status
      filter(n_all > 1) %>%
      
      mutate(status = case_when(
        n_mat_methylated == n_all ~ 'Maternal',
        n_bi_methylated > 0  & n_mat_methylated > 0~ 'Polymorphic',
      ))

    sandel_confirmed 

    ## # A tibble: 29 x 7
    ##    Gene  description n_mat_methylated n_pat_methylated n_bi_methylated n_all
    ##    <chr> <chr>                  <dbl>            <dbl>           <dbl> <dbl>
    ##  1 FRMD3 2 maternal…                2               NA               3     5
    ##  2 KCNQ1 8 maternal…                8               NA              NA     8
    ##  3 TMEM… 2 maternal…                2               NA              NA     2
    ##  4 TET3  5 maternal…                5               NA               3     8
    ##  5 SPHK… 4 maternal…                4               NA               5     9
    ##  6 ZNF3… 2 maternal…                2               NA              NA     2
    ##  7 C3OR… 2 maternal…                2               NA               1     3
    ##  8 FGF12 2 maternal…                2               NA               3     5
    ##  9 PDE6B 3 maternal…                3               NA              NA     3
    ## 10 STX1… 3 maternal…                3               NA              NA     3
    ## # … with 19 more rows, and 1 more variable: status <chr>

    sandel_confirmed <- sandel_confirmed %>% 
      mutate(methylated_allele = if_else(n_bi_methylated == n_all, 'bi', 'm', missing = 'm'),
             percent_biallelic = n_bi_methylated*100 / n_all) %>%
      select(Gene, methylated_allele, percent_biallelic, n_all) %>%
      dplyr::rename(n_used_to_verify_allele_specific_methylation = n_all) %>%
      
      # remove 100% biallelic
      filter(methylated_allele != 'bi')
    sandel_confirmed

    ## # A tibble: 22 x 4
    ##    Gene     methylated_allele percent_biallel… n_used_to_verify_allele_specific…
    ##    <chr>    <chr>                        <dbl>                             <dbl>
    ##  1 FRMD3    m                             60                                   5
    ##  2 KCNQ1    m                             NA                                   8
    ##  3 TMEM247  m                             NA                                   2
    ##  4 TET3     m                             37.5                                 8
    ##  5 SPHKAP   m                             55.6                                 9
    ##  6 ZNF385D  m                             NA                                   2
    ##  7 C3ORF62  m                             33.3                                 3
    ##  8 FGF12    m                             60                                   5
    ##  9 PDE6B    m                             NA                                   3
    ## 10 STX18-A… m                             NA                                   3
    ## # … with 12 more rows

column `n_used_to_verify_allele_specific_methylation` is the \# of
samples used to verify allele specific methylation. column
`percent_biallelic` is the percentage of samples that show biallelic
methylation. If &gt; 0, indicates polymorphic.

    # Which in s3 are not found in table s2
    sandel_confirmed$Gene %in% sandel_candidates$Gene

    ##  [1]  TRUE  TRUE  TRUE FALSE  TRUE  TRUE FALSE  TRUE FALSE  TRUE FALSE  TRUE
    ## [13]  TRUE FALSE FALSE  TRUE FALSE FALSE  TRUE FALSE FALSE  TRUE

    sandel_confirmed %>%
      filter(!Gene %in% sandel_candidates$Gene)

    ## # A tibble: 10 x 4
    ##    Gene     methylated_allele percent_biallel… n_used_to_verify_allele_specific…
    ##    <chr>    <chr>                        <dbl>                             <dbl>
    ##  1 TET3     m                             37.5                                 8
    ##  2 C3ORF62  m                             33.3                                 3
    ##  3 PDE6B    m                             NA                                   3
    ##  4 GPR78    m                             50                                   2
    ##  5 DENND3   m                             44.4                                 9
    ##  6 CACNA1C  m                             NA                                   2
    ##  7 PAPLN-A… m                             NA                                   2
    ##  8 Chr 18   m                             NA                                   3
    ##  9 ACTL10   m                             NA                                   2
    ## 10 TPTEP1   m                             NA                                   5

    # join based on "Gene"
    sandel <- sandel_candidates %>%
      fuzzyjoin::regex_right_join(sandel_confirmed, by = 'Gene') %>%
      select(contains('.x'), contains('.y'), everything())  %>%
      
      # remove s3 entries with no match in s2
      filter(!is.na(Gene.x)) %>%
      dplyr::rename(gene = Gene.y) %>%
      select(-Gene.x) %>%
      
      janitor::clean_names() %>%
      rename(chr = chromosome, cpg_island = cp_g_island, cpg_n =cp_gn)

    sandel

    ## # A tibble: 16 x 11
    ##    gene  chr    start    end width cpg_island cpg_n gene_pos methylated_alle…
    ##    <chr> <chr>  <dbl>  <dbl> <dbl> <chr>      <dbl> <chr>    <chr>           
    ##  1 FRMD3 chr9  8.61e7 8.61e7  1278 no            31 Gene_bo… m               
    ##  2 FRMD3 chr9  8.62e7 8.62e7  2911 yes          169 Promoter m               
    ##  3 KCNQ1 chr11 2.49e6 2.49e6  1722 no            53 Gene_bo… m               
    ##  4 TMEM… chr2  4.67e7 4.67e7  1491 no            52 Promoter m               
    ##  5 TET3  chr2  7.43e7 7.43e7  2712 yes           77 intrege… m               
    ##  6 SPHK… chr2  2.29e8 2.29e8  3526 yes          126 Promoter m               
    ##  7 ZNF3… chr3  2.18e7 2.18e7  1665 no            25 Promoter m               
    ##  8 FGF12 chr3  1.92e8 1.92e8  2824 yes          149 Promoter m               
    ##  9 STX1… chr4  4.58e6 4.58e6  1692 yes           86 Gene_bo… m               
    ## 10 SFRP2 chr4  1.55e8 1.55e8  6021 yes          312 Promoter m               
    ## 11 R3HC… chr8  2.31e7 2.31e7  1322 yes           59 Promoter m               
    ## 12 DENN… chr8  1.42e8 1.42e8  1288 yes           97 intrege… m               
    ## 13 CACN… chr12 2.80e6 2.80e6   358 yes           33 Promoter m               
    ## 14 FGF14 chr13 1.03e8 1.03e8  1856 yes          161 Promoter m               
    ## 15 CACN… chr19 1.36e7 1.36e7  3305 yes          140 Promoter m               
    ## 16 CACN… chr22 4.01e7 4.01e7  3728 yes          198 Gene_bo… m               
    ## # … with 2 more variables: percent_biallelic <dbl>,
    ## #   n_used_to_verify_allele_specific_methylation <dbl>

    sandel %>%
      write_tsv(here::here('processed', 'sanchez_delgado_final.tsv'))

Hamada 2016
-----------

**Genome Build:** hg19 **Tissues:** Placenta, Cytotrophoblasts, germ
cells, blood, stromal from CV **Technology:** WGBSEQ, RNAseq
**Approach:** Determine allele-specific meth/rna via sequencing in Cyto.
Supported by germ cell DMR analysis. **Findings:** 440 confirmed
maternally methylated primary imprinted DMRs (mDMRs). Out of 101 mDMRs
that Courtney identified, 31 were tested in this study, and 25 of those
were confirmed. An advantage of this study over previous is they used
entirely wgBSEQ -&gt; they suggest there are around 1800 imprinted DMRs
based on their results.

Limitation is discovery was an N of 1. (Candidates were followed up with
more samples in a targeted approach). Suggests a proportion of DMRs are
not placental-specific.

    hamada_s3 <- read_excel(here::here('raw files', 'Hamada 2016 s3.xlsx'))
    hamada_s4 <- read_excel(here::here('raw files', 'Hamada 2016 s4.xlsx'))

    # remove data columns
    hamada_s3 <- hamada_s3 %>%
      select(Classification:`Ranking according to mean [M-P]`) %>%
      mutate(Trimester = 'First')
    hamada_s4 <- hamada_s4 %>%
      select(Classification:`Ranking according to mean [M-P]`) %>%
      mutate(Trimester = 'Second/Term')

    hamada <- bind_rows(hamada_s3, hamada_s4) %>%
      janitor::clean_names()

Filter to significant DMRs, as defined in source:

> We defined candidate mDMRs or pDMRs showing &gt;30% or &lt;-30% mean
> \[M - P\] levels and statistically significant allelic methylation
> differences (BH-corrected p &lt; 0.05, Student’s t test) as confirmed
> gDMRs.

    hamada %>%
      count(trimester, classification) # 797 and 449 for 1st vs term/2nd

    ## # A tibble: 4 x 3
    ##   trimester   classification     n
    ##   <chr>       <chr>          <int>
    ## 1 First       Candidate mDMR   797
    ## 2 First       Candidate pDMR    97
    ## 3 Second/Term Candidate mDMR   449
    ## 4 Second/Term Candidate pDMR    43

    hamada %>%
      filter(corrected_p_value < 0.05, abs(mean_m_p_percent) >= 30) %>%
      count(trimester, classification)

    ## # A tibble: 4 x 3
    ##   trimester   classification     n
    ##   <chr>       <chr>          <int>
    ## 1 First       Candidate mDMR   384
    ## 2 First       Candidate pDMR    13
    ## 3 Second/Term Candidate mDMR   151
    ## 4 Second/Term Candidate pDMR     1

    hamada <- hamada %>%
      filter(corrected_p_value < 0.05, abs(mean_m_p_percent) >= 30) 

    hamada %>%
      write_tsv(here::here('processed', 'hamada_final.tsv'))

Terminology
===========

<table>
<colgroup>
<col style="width: 50%" />
<col style="width: 50%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Term</th>
<th style="text-align: left;">Definition</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">DMR</td>
<td style="text-align: left;">Sometimes imprinted regions characterized by allele-specific methylation are referred to as differentially methylated regions (DMR), where differential methylation refers to the opposite methylation states of the parental alleles</td>
</tr>
</tbody>
</table>

Merge
=====
