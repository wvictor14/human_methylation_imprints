This is a combined resource of 5 publications that looked at human
imprinting:

<table>
<thead>
<tr class="header">
<th>Publication</th>
<th>Tissue-specificity</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>Court et al. 2014</td>
<td>Placental + Other</td>
</tr>
<tr class="even">
<td>Hanna et al. 2016</td>
<td>Placental + Other</td>
</tr>
<tr class="odd">
<td>Sanchez-Delgado et al. 2016</td>
<td>Placental</td>
</tr>
<tr class="even">
<td>Hamada et al. 2016</td>
<td>Placental</td>
</tr>
<tr class="odd">
<td>Zink et al. 2018</td>
<td>Other</td>
</tr>
</tbody>
</table>

These regions can be found in @ `processed/all_imprinted_dmrs.tsv`

Here is a view of the first 6 rows:

<div class="kable-table">

<table>
<thead>
<tr>
<th style="text-align:right;">
chr
</th>
<th style="text-align:right;">
start
</th>
<th style="text-align:right;">
end
</th>
<th style="text-align:left;">
methylated\_allele
</th>
<th style="text-align:left;">
tissue\_specificity
</th>
<th style="text-align:left;">
associated\_gene
</th>
<th style="text-align:left;">
court
</th>
<th style="text-align:left;">
hanna
</th>
<th style="text-align:left;">
sanchez\_delgado
</th>
<th style="text-align:left;">
zink
</th>
<th style="text-align:left;">
hamada
</th>
<th style="text-align:left;">
note
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
6684860
</td>
<td style="text-align:right;">
6685996
</td>
<td style="text-align:left;">
M
</td>
<td style="text-align:left;">
placental-specific
</td>
<td style="text-align:left;">
THAP3
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
TRUE
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
</tr>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
7827139
</td>
<td style="text-align:right;">
7827709
</td>
<td style="text-align:left;">
M
</td>
<td style="text-align:left;">
other
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
TRUE
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
</tr>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
19614429
</td>
<td style="text-align:right;">
19615702
</td>
<td style="text-align:left;">
M
</td>
<td style="text-align:left;">
placental-specific
</td>
<td style="text-align:left;">
AKR7A3
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
TRUE
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
</tr>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
36184400
</td>
<td style="text-align:right;">
36184863
</td>
<td style="text-align:left;">
M
</td>
<td style="text-align:left;">
placental-specific
</td>
<td style="text-align:left;">
C1orf216
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
TRUE
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
</tr>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
38200920
</td>
<td style="text-align:right;">
38201123
</td>
<td style="text-align:left;">
M
</td>
<td style="text-align:left;">
other
</td>
<td style="text-align:left;">
EPHA10
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
TRUE
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
</tr>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
39559602
</td>
<td style="text-align:right;">
39559980
</td>
<td style="text-align:left;">
M
</td>
<td style="text-align:left;">
other
</td>
<td style="text-align:left;">
PPIEL
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
TRUE
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
</tr>
</tbody>
</table>
</div>

There are in total 418 imprinted regions, which can be divided into
those that are specific to the placenta, and those that have been
observed in other tissues.

<div class="kable-table">

<table>
<thead>
<tr>
<th style="text-align:left;">
tissue\_specificity
</th>
<th style="text-align:right;">
n
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
other
</td>
<td style="text-align:right;">
307
</td>
</tr>
<tr>
<td style="text-align:left;">
placental-specific
</td>
<td style="text-align:right;">
111
</td>
</tr>
</tbody>
</table>
</div>
