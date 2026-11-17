#!/usr/bin/env python3
"""patch_ground_truth_9_18b.py

I take TICK_PROJECT_GROUND_TRUTH_9_18.docx and write the 9_18 rev b version.

WHY A SECOND REVISION ON THE SAME DAY
-------------------------------------
The 9_18 file was written at midday. The afternoon overturned two of its own
new claims and one of its conclusions, so leaving it as the current ground
truth would leave a retracted result standing as the record.

WHAT CHANGED AFTER 9_18 WAS WRITTEN
-----------------------------------
    analysis unit          602 rows are 413 ANIMALS; 291 rows share an animal
    pooled co-occurrence   RETRACTED as a pathogen finding; it is one plate
    pool-size conditioning  invalid, wrong unit and NumTicks is not pool size
    burden conditioning    the burden null is 1.00; burden induces nothing
    detection vs pool size  flat, 1 tick detects as often as 40. OPEN.
    PoolTestR on 16S       invalid, whichever way the flat curve is explained

USAGE
-----
    python3 patch_ground_truth_9_18b.py IN.docx OUT.docx
"""
import os
import re
import shutil
import subprocess
import sys
import tempfile

APTOS = ('<w:rFonts w:ascii="Aptos" w:cs="Aptos" w:eastAsia="Aptos" '
         'w:hAnsi="Aptos"/>')

RPR_BODY = ('<w:rPr>' + APTOS + '<w:b w:val="false"/><w:bCs w:val="false"/>'
            '{ital}<w:color w:val="000000"/><w:sz w:val="22"/>'
            '<w:szCs w:val="22"/></w:rPr>')
RPR_CELL = ('<w:rPr>' + APTOS + '<w:b w:val="false"/><w:bCs w:val="false"/>'
            '{ital}<w:color w:val="000000"/><w:sz w:val="20"/>'
            '<w:szCs w:val="20"/></w:rPr>')
RPR_HEAD_CELL = ('<w:rPr>' + APTOS + '<w:b/><w:bCs/><w:color w:val="FFFFFF"/>'
                 '<w:sz w:val="20"/><w:szCs w:val="20"/></w:rPr>')
RPR_H2 = ('<w:rPr>' + APTOS + '<w:b/><w:bCs/><w:color w:val="6B3A4E"/>'
          '<w:sz w:val="30"/><w:szCs w:val="30"/></w:rPr>')
RPR_H3 = ('<w:rPr>' + APTOS + '<w:b/><w:bCs/>{ital}<w:color w:val="7D5A73"/>'
          '<w:sz w:val="25"/><w:szCs w:val="25"/></w:rPr>')
RPR_LIST = ('<w:rPr>' + APTOS + '{ital}<w:sz w:val="22"/>'
            '<w:szCs w:val="22"/></w:rPr>')

PPR_BODY = '<w:pPr><w:spacing w:after="120" w:line="276" w:lineRule="auto"/></w:pPr>'
PPR_H2 = ('<w:pPr><w:pStyle w:val="Heading2"/><w:keepNext/><w:keepLines/>'
          '<w:spacing w:after="120" w:before="320"/></w:pPr>')
PPR_H3 = ('<w:pPr><w:pStyle w:val="Heading3"/><w:keepNext/><w:keepLines/>'
          '<w:spacing w:after="100" w:before="240"/></w:pPr>')
PPR_LIST = ('<w:pPr><w:pStyle w:val="ListParagraph"/><w:numPr>'
            '<w:ilvl w:val="0"/><w:numId w:val="2"/></w:numPr>'
            '<w:spacing w:after="80"/></w:pPr>')

TCPR = ('<w:tcPr><w:tcW w:type="dxa" w:w="0"/><w:shd w:fill="{fill}" '
        'w:val="clear"/><w:tcMar><w:top w:type="dxa" w:w="80"/>'
        '<w:left w:type="dxa" w:w="100"/><w:bottom w:type="dxa" w:w="80"/>'
        '<w:right w:type="dxa" w:w="100"/></w:tcMar></w:tcPr>')


def escape(text):
    return (text.replace('&', '&amp;').replace('<', '&lt;')
                .replace('>', '&gt;'))


def runs(text, rpr_template):
    """I turn text carrying *asterisk* italic spans into w:r runs."""
    out = []
    for index, piece in enumerate(text.split('*')):
        if piece == '':
            continue
        ital = ('<w:i/><w:iCs/>' if index % 2
                else '<w:i w:val="false"/><w:iCs w:val="false"/>')
        out.append('<w:r>' + rpr_template.format(ital=ital)
                   + '<w:t xml:space="preserve">' + escape(piece)
                   + '</w:t></w:r>')
    return ''.join(out)


def body_paragraph(text):
    return '<w:p>' + PPR_BODY + runs(text, RPR_BODY) + '</w:p>'


def heading2(text):
    return ('<w:p>' + PPR_H2 + '<w:r>' + RPR_H2
            + '<w:t xml:space="preserve">' + escape(text) + '</w:t></w:r></w:p>')


def heading3(text):
    """I build a Heading 3, honouring *asterisk* italics for species names."""
    return '<w:p>' + PPR_H3 + runs(text, RPR_H3) + '</w:p>'


def list_item(text):
    return '<w:p>' + PPR_LIST + runs(text, RPR_LIST) + '</w:p>'


def header_row(cells):
    out = ''
    for text in cells:
        out += ('<w:tc>' + TCPR.format(fill='4A2233') + '<w:p><w:r>'
                + RPR_HEAD_CELL + '<w:t xml:space="preserve">' + escape(text)
                + '</w:t></w:r></w:p></w:tc>')
    return '<w:tr><w:trPr><w:tblHeader/></w:trPr>' + out + '</w:tr>'


def data_row(cells, fill):
    out = ''
    for text in cells:
        out += ('<w:tc>' + TCPR.format(fill=fill) + '<w:p>'
                + runs(text, RPR_CELL) + '</w:p></w:tc>')
    return '<w:tr>' + out + '</w:tr>'


def replace_once(xml, old, new, what):
    count = xml.count(old)
    if count != 1:
        sys.exit('expected exactly 1 occurrence of ' + what + ', found '
                 + str(count) + '. The source document is not the one this '
                 'script was written against; do not trust the output.')
    return xml.replace(old, new, 1)


def replace_text(xml, old_text, new_text, what):
    """I replace one whole run's visible text."""
    return replace_once(xml, escape(old_text), escape(new_text), what)



# ---------------------------------------------------------------- new content

SECTION_1_ADDITIONS = [
    heading3('1.21  The unit of analysis is the ANIMAL, and the analysis table '
             'was keyed on the library.'),
    body_paragraph(
        'The 602 rows of analysis_table_16S.tsv are 413 animals. 102 animals '
        'hold more than one row, accounting for 291 of 602 rows (48.3%). The '
        'Barcode column and the library lookup produce the same 413 groups, so '
        'this is not a resolver artifact.'),
    body_paragraph(
        '602 was never chosen as a sample size. build_analysis_table.py emits '
        'one row per 16S library because nothing instructed it otherwise, and '
        '602 became "n" by being the number of rows in a file. The design is a '
        'pooled design in which the POOL IS THE ANIMAL: 413 pools, pool size '
        'NumTicks, and the 602 rows are repeat measurements of those 413 '
        'pools.'),
    body_paragraph(
        'NumTicks is the animal\'s tick count copied onto each of its rows. '
        'Summed over rows it is 4,730; summed over animals it is 3,209, '
        'against 3,939 in the crosswalk across all 563 specimens. Every one of '
        'the 102 repeated animals carries a single value across all its rows.'),
    body_paragraph(
        'The intent was recorded nowhere a script could see it. This document '
        'had no unit-of-analysis section before rev b, which is how the drift '
        'went unnoticed for the life of the project.'),

    heading3('1.22  Collapsing with any-positive is BIASED and must not be '
             'used.'),
    body_paragraph(
        'An animal counted positive if any of its rows is positive gets one '
        'detection opportunity per library. Positivity by rows per animal: '
        '*Anaplasma* 38.3% at one row, 59.0% at two, 64.7% at three, 83.3% at '
        'four or more. *Borreliella* 37.0%, 52.5%, 52.9%, 66.7%.'),
    body_paragraph(
        'That collapse raises *Anaplasma* from 38.9% (234/602) to 45.0% '
        '(186/413) and takes *Anaplasma* with *Borrelia* from 3.08 to 8.00. '
        'Both are artifacts.'),
    body_paragraph(
        'The defensible animal-level estimates come from drawing one row per '
        'animal at random, 2,000 times, seed 20260918: *Anaplasma* 38.5% (36.8 '
        'to 40.0), *Borreliella* 35.6% (33.9 to 37.3), *Borrelia* 3.6% (2.9 to '
        '4.4), *Bartonella* 13.1% (12.1 to 14.3). On the 311 animals holding '
        'exactly one row, *Anaplasma* is 38.3% (119/311, 95% CI 33.0 to 43.8) '
        'and *Borreliella* 37.0% (115/311, 31.8 to 42.5).'),
    body_paragraph(
        'Table 1 percentages move very little. The denominator moves from 602 '
        'to 413 and the intervals widen accordingly.'),

    heading3('1.23  The pooled co-occurrence is confined to ONE sequencing '
             'plate. RETRACTED as a pathogen finding.'),
    body_paragraph(
        '*Anaplasma* with *Borreliella* survives the change of unit: 2.03 per '
        'row, 2.00 on the 311 single-library animals (p = 0.005), 1.96 '
        'resampled across 413 (1.68 to 2.29), significant in 2,000 of 2,000 '
        'draws. It is not the same association in every plate.'),
    body_paragraph(
        'Breslow-Day at animal level rejects homogeneity with median p = '
        '0.024, in 85% of 1,000 draws, against a measured 7% false-positive '
        'rate. 241209_P3-5 runs at 6.80 against roughly 1.5 elsewhere, and '
        'dropping it takes the crude from 1.89 to 1.50.'),
    body_paragraph(
        'Against a null built on the measured burden curve (1.25): all 413 '
        'animals 1.97 against a null of 1.00, p = 0.001 · excluding P3-5, 1.47 '
        'against 1.01, p = 0.053 on 328 animals · P3-5 alone, 6.92 against '
        '1.05, p < 0.001 on 85 animals. *Borreliella* with *Bartonella* '
        'behaves the same way: 2.06 overall, 1.44 outside P3-5 at p = 0.127, '
        '4.94 inside at p = 0.018.'),
    body_paragraph(
        'Outside that one plate NOTHING reaches significance across 328 '
        'animals. The pooled co-occurrence is a property of 241209_P3-5 and is '
        'not reportable as a finding about pathogens.'),
    body_paragraph(
        '*Anaplasma* with *Borrelia* should be reported as a count, not an '
        'odds ratio. There are 23 *Borrelia* positives in 602 rows and the '
        'resampled interval runs from 3.49 to 18.80.'),

    heading3('1.24  Tick burden induces NO association in this dataset.'),
    body_paragraph(
        'Conditioning the animal-level co-occurrence on the animal\'s own tick '
        'count, using a null built from the observed positivity at each tick '
        'band, gives a null odds ratio of 1.00. The argument that an animal '
        'carrying more ticks has more chances to show both organisms is '
        'empirically empty here, and every version of it from 2026-09-18, in '
        'both the pool-size and the burden form, was wrong.'),
    body_paragraph(
        'Quartile stratification does not show this. On simulated animals '
        'whose organisms are independent tick by tick, a crude 3.32 left 1.50 '
        'after quartile adjustment, because burden varies four-fold inside the '
        'top quartile. Quartile-adjusted odds ratios from 2026-09-18 should '
        'not be quoted.'),

    heading3('1.25  Detection does NOT increase with pool size. OPEN.'),
    body_paragraph(
        'On the 311 single-library animals, where the pool size is the '
        'animal\'s tick count and is not in dispute: *Borreliella* 38.9% '
        '(28/72) at one tick, 43.9% (29/66) at two, 36.0% (9/25) at 27 to 62. '
        '*Bartonella* 11.1% (8/72) at one tick, 12.0% (3/25) at 27 to 62. '
        '*Anaplasma* 26.4% (19/72) at one tick, 44.0% (11/25) at 27 to 62.'),
    body_paragraph(
        'A binomial pool model fitted to the same data predicts 87.4%, 44.9% '
        'and 91.0% respectively at the top band. A pool of forty ticks detects '
        'at the same rate as a single tick.'),
    body_paragraph(
        'TWO READINGS, NOT YET DISTINGUISHED. Infection may be a property of '
        'the ANIMAL, so that ticks off one deer share its exposure and a pool '
        'of forty reports what the deer had; that is a reportable finding. Or '
        'NumTicks does not describe what went into the tube, which is a '
        'records problem.'),
    body_paragraph(
        'The within-animal agreement test cannot separate them from the files '
        'on disk. All three single-tick tables hold one tick per animal, 93 '
        'ticks from 93 animals. The only table with repeats gives 13 animals '
        'and 25 pairs, where one of six tests reached p = 0.047, which is what '
        'chance produces across six tests.'),
    body_paragraph(
        'WHAT DECIDES IT: microbiome richness against pool size on the '
        'single-library animals, conditioned on read depth, since a library '
        'built from forty ticks should carry more ASVs than one built from '
        'one. Or per-tick qPCR results for several ticks from the same animal.'),

    heading3('1.26  PoolTestR cannot be applied to the 16S data.'),
    body_paragraph(
        'Its estimator infers per-tick prevalence from pool size on the '
        'assumption that ticks within a pool are independent and that '
        'detection does not depend on pool size. 1.25 shows pool size carries '
        'no information about detection in this assay. 38.9% positive at one '
        'tick implies a per-tick *Borreliella* prevalence near 0.39; 36.0% at '
        '27 to 62 ticks implies about 0.014. Those differ more than twenty-fold '
        'and the pooled figure is the wrong one.'),
    body_paragraph(
        'This holds whichever reading of 1.25 is correct, and it applies to '
        'the 16S data. Applying the estimator to per-tick qPCR results is a '
        'separate question and is not affected.'),
    body_paragraph(
        'The per-tick column of Table 1 in results/validation/'
        'pooltestr_validation_numbers.md, sent to Joe and Lawrence in August, '
        'rests on this estimator and needs revisiting.'),

    heading3('1.27  RunBatch is unreliable; the run token is not.'),
    body_paragraph(
        'RunBatch is blank on 43 of 602 rows and contradicts the sample name '
        'on 4. The 143 rows whose sample id ends 241209_P3-5 carry: '
        'AMP-120924-USGS-TICKS-P3-5 on 96, nothing on 43, and '
        'AMP-072924-DN-USGS-Tick-2024-07-16-Batch-02-I8 on 4. Those 4 are '
        'ALHT-Deer-64, ALHT-Deer-92, ALHT-Deer-93-1 and ALHT-Deer-93-2, and '
        'the run sheet needs to say which label is right.'),
    body_paragraph(
        'The run token at the end of every sample id is present on all 602 '
        'rows. Batch stratification should use the token.'),
    body_paragraph(
        '241209_P3-5 is 143 rows from 85 animals, 1.68 rows per animal against '
        'about 1.4 elsewhere, hosts 75 deer, 8 moose and 2 bear, tick counts '
        'from 1 to 54 with a median of 5. It is the only plate mixing three '
        'host species.'),
]

RETRACTION = [
    body_paragraph(
        'RETRACTED 2026-09-18 rev b. Everything in this block was computed on '
        '602 rows treated as 602 independent pools. They are 413 animals '
        '(1.21). The pool-size conditioning is invalid twice over: wrong unit, '
        'and NumTicks is the animal\'s total repeated rather than the size of '
        'any tube. Conditioning on burden at the animal level gives a null '
        'odds ratio of 1.00, so burden induces nothing (1.24), and the '
        'association itself is confined to one plate (1.23). The numbers below '
        'are kept as the record of what was believed; do not quote them.'),
]

CLOSED_ADDITIONS = [
    list_item('The pooled co-occurrence question is CLOSED, against the '
              'expected answer. It is not a pathogen finding. It is one '
              'sequencing plate, 241209_P3-5 (1.23).'),
    list_item('Batch conditioning is done, at the animal level, and the '
              'batches do NOT share one odds ratio: Breslow-Day median p = '
              '0.024, rejecting in 85% of 1,000 draws.'),
    list_item('The analysis unit is settled: 413 animals, not 602 rows '
              '(1.21). Animal-level prevalence is by one-row-per-animal '
              'resampling, never by any-positive (1.22).'),
    list_item('Tick burden is tested and induces nothing (1.24). The pool-size '
              'argument that shaped a day of analysis is empirically empty in '
              'this data.'),
    list_item('Five read-only checkers were written and validated against '
              'synthetic data with the answer known BEFORE being run on the '
              'real table: check_specimen_replication.py, '
              'check_animal_level_unit.py, check_batch_at_animal_level.py, '
              'check_burden_empirical_null.py and '
              'check_within_animal_clustering.py.'),
]

OPEN_ADDITIONS = [
    list_item('1.25 IS THE LIVE QUESTION. Detection does not increase with '
              'pool size, and the two readings (animal-level infection versus '
              'a records fault in NumTicks) lead opposite ways. DECIDE IT with '
              'microbiome richness against pool size on the 311 single-library '
              'animals, conditioned on read depth. Richness that climbs with '
              'NumTicks means the pools are real and the flat curve is '
              'biology; richness that is flat too means the tubes do not hold '
              'what the metadata says.'),
    list_item('GET PER-TICK qPCR RESULTS. Several ticks from the same animal, '
              'measured separately, settle 1.25 directly and nothing else on '
              'disk can. The three single-tick tables are one tick per animal '
              'by design.'),
    list_item('ALIQUOTS OR DIVIDED, for the 102 multi-library animals. On the '
              'measured curve the two predictions land within a few points of '
              'each other, so the data cannot say. It is a protocol question: '
              'when an animal has N ticks and L subsample libraries, were the '
              'ticks divided among the libraries or was one pooled extract '
              'split? EQ02441291 is the case to ask Kellen about: 11 ticks, 10 '
              'subsample libraries.'),
    list_item('RunBatch is blank on 43 rows and contradicts the sample name on '
              '4 (1.27). Resolve ALHT-Deer-64, ALHT-Deer-92, ALHT-Deer-93-1 '
              'and ALHT-Deer-93-2 against the run sheet.'),
    list_item('WHAT IS 241209_P3-5? It carries the whole co-occurrence signal, '
              'it is the only plate mixing three host species, and it holds '
              'the subsampled animals. Until it is understood the plate '
              'difference cannot be attributed to anything, and batch cannot '
              'be separated from animal in this design anyway.'),
    list_item('results/cooccurrence/single_tick_cooccurrence.tsv has 124 rows '
              'against the current cohort of 93. Check whether it is stale '
              'before anything reads it.'),
    list_item('PoolTestR on 16S is retired (1.26). The August memo to Joe and '
              'Lawrence carries a per-tick column that rests on it.'),
]

PATTERN_ADDITIONS = [
    body_paragraph(
        'A SECOND SHAPE, added 2026-09-18 rev b. The first eight bugs were one '
        'key that looked unique and was not. These are a number computed on '
        'the wrong unit, or checked against an assumed model instead of the '
        'data sitting in the same output.'),
    list_item('An analysis keyed on the library when the design\'s unit is the '
              'animal, so 602 rows became an n without anyone choosing it.'),
    list_item('A collapse rule, any-positive, that created the bias it was '
              'meant to remove: more libraries means more chances to be '
              'called positive.'),
    list_item('A confounder conditioned on by quartiles too coarse to remove '
              'it, leaving 1.50 of a 3.32 that was entirely confounding, and '
              'the residue read as survival.'),
    list_item('A null built on an assumed model and reported as p = 1.000, '
              'while the same script\'s own output two sections earlier showed '
              'the model failing.'),
    list_item('A summary statistic quoted across strata that do not share a '
              'common value, where Mantel-Haenszel averages disagreeing strata '
              'into a number that describes none of them.'),
    body_paragraph(
        'The counter-practice that worked, and should be kept: measure a '
        'test\'s POWER before reading its result. Doing so stopped a clean '
        'batch result being reported from a test with 28% power, and it is why '
        'the one-plate-against-the-rest test exists at all.'),
]


def main():
    if len(sys.argv) != 3:
        sys.exit('usage: patch_ground_truth_9_18b.py IN.docx OUT.docx')
    source, target = sys.argv[1], sys.argv[2]
    if not os.path.isfile(source):
        sys.exit('input not found: ' + source)

    workdir = tempfile.mkdtemp(prefix='gt18b_')
    try:
        if subprocess.run(['unzip', '-q', '-o', source,
                           '-d', workdir]).returncode != 0:
            sys.exit('could not unpack ' + source)
        doc_path = os.path.join(workdir, 'word', 'document.xml')
        if not os.path.isfile(doc_path):
            sys.exit('no word/document.xml in ' + source)
        with open(doc_path, encoding='utf-8') as handle:
            xml = handle.read()

        if 'rev b' in xml:
            sys.exit('this document already says "rev b". Nothing changed.')

        # ---- revision marks -------------------------------------------------
        xml = replace_once(xml, '>2026-09-18<', '>2026-09-18 rev b<',
                           'the date line')
        xml = replace_text(
            xml,
            'TICK PROJECT GROUND TRUTH  |  2026-09-18  |  supersedes sections '
            '1, 2, 3, 9, 12 of the 8/6 file',
            'TICK PROJECT GROUND TRUTH  |  2026-09-18 rev b  |  supersedes '
            'sections 1, 2, 3, 9, 12 of the 8/6 file, and the pooled '
            'co-occurrence block of the 9/18 midday file',
            'the closing rule line')

        # ---- section 1 additions, after 1.20 --------------------------------
        anchor = escape('Nobody has used the pooled positive counts the sheet '
                        'contains.')
        hit = xml.find(anchor)
        if hit == -1:
            sys.exit('could not find the end of section 1.20')
        at = xml.find('</w:p>', hit) + len('</w:p>')
        xml = xml[:at] + ''.join(SECTION_1_ADDITIONS) + xml[at:]

        # ---- retract the pooled co-occurrence block -------------------------
        xml = replace_text(
            xml,
            'Pooled co-occurrence, n = 602, and what it survives',
            'Pooled co-occurrence, n = 602: RETRACTED, see 1.21 to 1.26',
            'the pooled co-occurrence heading')
        anchor = escape('Pooled co-occurrence, n = 602: RETRACTED, see '
                        '1.21 to 1.26')
        hit = xml.find(anchor)
        if hit == -1:
            sys.exit('could not find the pooled co-occurrence heading')
        at = xml.find('</w:p>', hit) + len('</w:p>')
        xml = xml[:at] + ''.join(RETRACTION) + xml[at:]

        # ---- the per-tick-null paragraph rested on it -----------------------
        xml = replace_text(
            xml,
            'DO NOT read that per-tick null as settling the question. At n = 93 '
            'it cannot. The pooled analysis on 602 samples gives the SAME pair '
            'an odds ratio of 2.03 at q = 0.0002, and that survives '
            'conditioning on pool size and read depth (see the pooled '
            'co-occurrence block below). OR 1.31 on 93 ticks and OR 2.05 on '
            '602 samples are the same direction with a seventh of the data; 8 '
            'co-detections against roughly 7 expected distinguishes nothing. '
            'An earlier draft of this document presented the two as '
            'contradicting each other. They do not.',
            'DO NOT read that per-tick null as settling the question. At n = 93 '
            'it cannot; 8 co-detections against roughly 7 expected '
            'distinguishes nothing. REVISED 2026-09-18 rev b: the pooled '
            'comparison that used to sit here is withdrawn. Those 602 rows are '
            '413 animals (1.21), the pooled association is confined to one '
            'plate (1.23), and outside that plate it does not reach '
            'significance across 328 animals. The per-tick analysis is now the '
            'only co-occurrence result, and it is underpowered rather than '
            'contradicted.',
            'the per-tick underpowered paragraph')

        # ---- section 12, Closed ---------------------------------------------
        xml = replace_text(xml, 'Closed 9 to 18 September 2026',
                           'Closed 9 to 18 September 2026 (rev b)',
                           'the Closed heading')
        anchor = escape('emitted by build_4marker_crosswalk.py, 9,409 keys, '
                        'zero collisions.')
        hit = xml.find(anchor)
        if hit == -1:
            sys.exit('could not find the end of the Closed list')
        at = xml.find('</w:p>', hit) + len('</w:p>')
        xml = xml[:at] + ''.join(CLOSED_ADDITIONS) + xml[at:]

        # ---- section 12, Open: the batch item is answered now ---------------
        stale = escape('CONDITION THE POOLED CO-OCCURRENCE ON BATCH.')
        hit = xml.find(stale)
        if hit == -1:
            sys.exit('could not find the batch item to close')
        item_start = xml.rfind('<w:p>', 0, hit)
        item_end = xml.find('</w:p>', hit) + len('</w:p>')
        xml = xml[:item_start] + xml[item_end:]

        anchor = escape('build_gt.js and build_postmortem.js DO NOT EXIST.')
        hit = xml.find(anchor)
        if hit == -1:
            sys.exit('could not find the head of the Still open list')
        at = xml.rfind('<w:p>', 0, hit)
        xml = xml[:at] + ''.join(OPEN_ADDITIONS) + xml[at:]

        # ---- section 14 -----------------------------------------------------
        anchor = escape('185 specimens, under a comment asserting the drops '
                        'were expected.')
        hit = xml.find(anchor)
        if hit == -1:
            sys.exit('could not find the end of the section 14 list')
        at = xml.find('</w:p>', hit) + len('</w:p>')
        xml = xml[:at] + ''.join(PATTERN_ADDITIONS) + xml[at:]

        with open(doc_path, 'w', encoding='utf-8') as handle:
            handle.write(xml)

        footer_path = os.path.join(workdir, 'word', 'footer1.xml')
        if os.path.isfile(footer_path):
            with open(footer_path, encoding='utf-8') as handle:
                footer = handle.read()
            footer = footer.replace('Ground truth 2026-09-18',
                                    'Ground truth 2026-09-18 rev b')
            with open(footer_path, 'w', encoding='utf-8') as handle:
                handle.write(footer)

        if os.path.exists(target):
            os.remove(target)
        rc = subprocess.run(['zip', '-q', '-r', '-X', os.path.abspath(target),
                             '.'], cwd=workdir).returncode
        if rc != 0:
            sys.exit('could not repack to ' + target)
        print('wrote ' + target)
        print()
        print('New: 1.21 unit is the animal · 1.22 any-positive is biased ·')
        print('1.23 co-occurrence RETRACTED, it is one plate · 1.24 burden')
        print('induces nothing · 1.25 flat detection curve, OPEN · 1.26')
        print('PoolTestR retired for 16S · 1.27 RunBatch unreliable.')
        print('Section 2 pooled block marked RETRACTED, section 12 Closed and')
        print('Open updated, section 14 second shape added.')
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


if __name__ == '__main__':
    main()
