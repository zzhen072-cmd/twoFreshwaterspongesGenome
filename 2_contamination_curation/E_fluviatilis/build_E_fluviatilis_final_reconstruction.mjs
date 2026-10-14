#!/usr/bin/env node

import { createHash } from 'node:crypto';
import { createReadStream, createWriteStream } from 'node:fs';
import { mkdir, readFile, rename, stat, writeFile } from 'node:fs/promises';
import { dirname, join, resolve } from 'node:path';
import { createInterface } from 'node:readline';
import { finished } from 'node:stream/promises';
import { fileURLToPath } from 'node:url';
import { createGunzip, createGzip } from 'node:zlib';

// The archived execution used a local absolute path.  The public copy accepts
// an explicit project root and otherwise uses the current working directory.
const ROOT = resolve(process.env.FRESHWATER_PROJECT_ROOT ?? process.cwd());
const OUT_DIR = join(ROOT, 'work/nj151201_final_reconstruction');
const SRC_DIR = join(ROOT, 'output/manuscript_validation_20260923/annotation/deposit');
const DEFAULT_DELETION_LIST = fileURLToPath(
  new URL('./E_fluviatilis_removed_21_contigs.txt', import.meta.url),
);

const SOURCES = {
  genome: join(ROOT, 'output/manuscript_validation_20260923/assemblies/S04.genome.fna.gz'),
  deletionList: resolve(
    process.env.E_FLUVIATILIS_REMOVAL_LIST ?? DEFAULT_DELETION_LIST,
  ),
  gff: join(SRC_DIR, 'S04.annotation.original_contig_ids.gff3.gz'),
  cdsAll: join(SRC_DIR, 'S04.CDS.all_isoforms.fna.gz'),
  cdsRep: join(SRC_DIR, 'S04.CDS.representative.fna.gz'),
  proteinAll: join(SRC_DIR, 'S04.proteins.all_isoforms.faa.gz'),
  proteinRep: join(SRC_DIR, 'S04.proteins.representative.faa.gz'),
  mapping: join(SRC_DIR, 'S04_contig_id_mapping.tsv'),
  functional: join(SRC_DIR, 'S04.functional_and_model_annotations.tsv.gz'),
};

const OUTPUTS = {
  genome: join(OUT_DIR, 'NJ151201.genome.fna.gz'),
  gff: join(OUT_DIR, 'NJ151201.annotation.original_contig_ids.gff3.gz'),
  cdsAll: join(OUT_DIR, 'NJ151201.CDS.all_isoforms.fna.gz'),
  cdsRep: join(OUT_DIR, 'NJ151201.CDS.representative.fna.gz'),
  proteinAll: join(OUT_DIR, 'NJ151201.proteins.all_isoforms.faa.gz'),
  proteinRep: join(OUT_DIR, 'NJ151201.proteins.representative.faa.gz'),
  mapping: join(OUT_DIR, 'NJ151201.contig_id_mapping.tsv.gz'),
  functional: join(OUT_DIR, 'NJ151201.functional_and_model_annotations.tsv.gz'),
};

const EXPECTED = {
  deletedContigs: 21,
  contigs: 2097,
  assemblyBp: 267354623,
  n50: 237588,
  n90: 62222,
  l50: 315,
  l90: 1151,
  genes: 29180,
  transcripts: 33684,
  representatives: 29180,
  functional: {
    any_upset_source: 18082,
    eggNOG: 13267,
    InterPro: 16541,
    SwissProt: 11644,
    GO_eggNOG: 8566,
    GO_raw_union: 14111,
    KEGG_KO_union: 8649,
    COG_categories: 12431,
  },
};

function assert(condition, message) {
  if (!condition) throw new Error(`VALIDATION FAILED: ${message}`);
}

async function* linesFromPlain(path) {
  const input = createReadStream(path);
  const rl = createInterface({ input, crlfDelay: Infinity });
  for await (const line of rl) yield line;
}

async function* linesFromGzip(path) {
  const input = createReadStream(path).pipe(createGunzip());
  const rl = createInterface({ input, crlfDelay: Infinity });
  for await (const line of rl) yield line;
}

async function withAtomicGzip(path, producer) {
  await mkdir(dirname(path), { recursive: true });
  const tmp = `${path}.tmp`;
  const gzip = createGzip({ level: 9, mtime: 0 });
  const output = createWriteStream(tmp);
  gzip.pipe(output);
  try {
    await producer(gzip);
    gzip.end();
    await finished(output);
    await rename(tmp, path);
  } catch (error) {
    gzip.destroy(error);
    output.destroy(error);
    throw error;
  }
}

async function writeTextAtomic(path, text) {
  await mkdir(dirname(path), { recursive: true });
  const tmp = `${path}.tmp`;
  await writeFile(tmp, text, 'utf8');
  await rename(tmp, path);
}

async function sha256File(path) {
  const hash = createHash('sha256');
  for await (const chunk of createReadStream(path)) hash.update(chunk);
  return hash.digest('hex');
}

function sequenceHash(sequenceLines) {
  const hash = createHash('sha256');
  for (const line of sequenceLines) hash.update(line.trim().toUpperCase());
  return hash.digest('hex');
}

function fastaId(header) {
  assert(header.startsWith('>'), `FASTA header does not start with >: ${header.slice(0, 80)}`);
  return header.slice(1).trim().split(/\s+/, 1)[0];
}

async function filterFasta(source, destination, keepId, collectSequenceStats = false) {
  const sourceLines = linesFromGzip(source);
  const retainedIds = new Set();
  const allIds = new Set();
  const removedIds = new Set();
  const lengths = [];
  const sequenceSha256 = new Map();
  let sourceRecords = 0;
  let retainedRecords = 0;
  let header = null;
  let sequenceLines = [];

  await withAtomicGzip(destination, async (writer) => {
    const flush = () => {
      if (header === null) return;
      const id = fastaId(header);
      assert(!allIds.has(id), `duplicate FASTA ID in ${source}: ${id}`);
      allIds.add(id);
      sourceRecords += 1;
      if (keepId(id)) {
        retainedIds.add(id);
        retainedRecords += 1;
        writer.write(`${header}\n`);
        for (const line of sequenceLines) writer.write(`${line}\n`);
        if (collectSequenceStats) {
          const length = sequenceLines.reduce((sum, line) => sum + line.trim().length, 0);
          lengths.push(length);
          sequenceSha256.set(id, sequenceHash(sequenceLines));
        }
      } else {
        removedIds.add(id);
      }
    };

    for await (const line of sourceLines) {
      if (line.startsWith('>')) {
        flush();
        header = line;
        sequenceLines = [];
      } else {
        assert(header !== null, `sequence before first FASTA header in ${source}`);
        sequenceLines.push(line);
      }
    }
    flush();
  });

  return {
    sourceRecords,
    retainedRecords,
    allIds,
    retainedIds,
    removedIds,
    lengths,
    sequenceSha256,
  };
}

function parseAttributes(text) {
  const result = new Map();
  for (const token of text.split(';')) {
    if (!token) continue;
    const index = token.indexOf('=');
    if (index < 0) continue;
    result.set(token.slice(0, index), token.slice(index + 1));
  }
  return result;
}

function nxStats(lengths, fraction) {
  const sorted = [...lengths].sort((a, b) => b - a);
  const total = sorted.reduce((sum, n) => sum + n, 0);
  const threshold = total * fraction;
  let cumulative = 0;
  for (let i = 0; i < sorted.length; i += 1) {
    cumulative += sorted[i];
    if (cumulative >= threshold) return { nx: sorted[i], lx: i + 1 };
  }
  throw new Error('cannot compute Nx for an empty assembly');
}

async function filterGff(source, destination, deletedContigs) {
  const genes = new Set();
  const transcripts = new Set();
  const transcriptToGene = new Map();
  const removedGenes = new Set();
  const removedTranscripts = new Set();
  const keptFeatureCounts = {};
  const removedSeqidsSeen = new Set();

  await withAtomicGzip(destination, async (writer) => {
    for await (const line of linesFromGzip(source)) {
      if (line === '' || line.startsWith('#')) {
        writer.write(`${line}\n`);
        continue;
      }
      const fields = line.split('\t');
      assert(fields.length === 9, `GFF3 row does not have 9 fields: ${line.slice(0, 160)}`);
      const [seqid, , type, , , , , , attributeText] = fields;
      const attrs = parseAttributes(attributeText);
      const removed = deletedContigs.has(seqid);
      if (removed) removedSeqidsSeen.add(seqid);

      if (type === 'gene') {
        const id = attrs.get('ID');
        assert(id, `gene row without ID on ${seqid}`);
        (removed ? removedGenes : genes).add(id);
      } else if (type === 'mRNA' || type === 'transcript') {
        const id = attrs.get('ID');
        const parent = attrs.get('Parent');
        assert(id && parent, `${type} row lacks ID or Parent on ${seqid}`);
        if (removed) {
          removedTranscripts.add(id);
        } else {
          transcripts.add(id);
          transcriptToGene.set(id, parent.split(',')[0]);
        }
      }

      if (!removed) {
        keptFeatureCounts[type] = (keptFeatureCounts[type] ?? 0) + 1;
        writer.write(`${line}\n`);
      }
    }
  });

  assert(removedSeqidsSeen.size === deletedContigs.size,
    `GFF3 contains ${removedSeqidsSeen.size}/${deletedContigs.size} deletion-list contigs`);
  for (const transcript of transcripts) {
    const gene = transcriptToGene.get(transcript);
    assert(genes.has(gene), `retained transcript ${transcript} points to absent gene ${gene}`);
  }

  return {
    genes,
    transcripts,
    transcriptToGene,
    removedGenes,
    removedTranscripts,
    keptFeatureCounts,
    removedSeqidsSeen,
  };
}

async function filterMapping(source, destination, deletedContigs, retainedGenome) {
  let retainedRows = 0;
  let removedRows = 0;
  const retainedOriginalIds = new Set();
  const removedOriginalIds = new Set();

  await withAtomicGzip(destination, async (writer) => {
    let header = null;
    let indexes = null;
    for await (const line of linesFromPlain(source)) {
      if (header === null) {
        header = line.split('\t');
        indexes = Object.fromEntries(header.map((name, index) => [name, index]));
        for (const required of ['short_id', 'original_id', 'length_bp', 'sequence_sha256']) {
          assert(indexes[required] !== undefined, `mapping lacks required column ${required}`);
        }
        writer.write(`${line}\n`);
        continue;
      }
      if (line === '') continue;
      const fields = line.split('\t');
      const originalId = fields[indexes.original_id];
      if (deletedContigs.has(originalId)) {
        removedRows += 1;
        removedOriginalIds.add(originalId);
        continue;
      }
      assert(retainedGenome.retainedIds.has(originalId), `mapping original_id absent from retained genome: ${originalId}`);
      assert(Number(fields[indexes.length_bp]) === retainedGenome.lengthById.get(originalId),
        `mapping length mismatch for ${originalId}`);
      assert(fields[indexes.sequence_sha256] === retainedGenome.sequenceSha256.get(originalId),
        `mapping sequence SHA256 mismatch for ${originalId}`);
      retainedRows += 1;
      retainedOriginalIds.add(originalId);
      writer.write(`${line}\n`);
    }
  });

  assert(removedRows === deletedContigs.size, `mapping removed ${removedRows}, expected ${deletedContigs.size}`);
  assert(removedOriginalIds.size === deletedContigs.size,
    `mapping contains ${removedOriginalIds.size}/${deletedContigs.size} deletion-list contigs`);
  assert(retainedOriginalIds.size === retainedGenome.retainedIds.size,
    `mapping retained ID count ${retainedOriginalIds.size} != genome ${retainedGenome.retainedIds.size}`);
  return { retainedRows, removedRows, retainedOriginalIds, removedOriginalIds };
}

function nonempty(value) {
  return value !== undefined && value !== '' && value !== '[]' && value !== 'false' && value !== 'False' && value !== '0';
}

async function filterFunctional(source, destination, retainedRepresentativeIds, deletedContigs) {
  let sourceRows = 0;
  let retainedRows = 0;
  let removedRows = 0;
  let header = null;
  let indexes = null;
  const retainedProteinIds = new Set();
  const counts = Object.fromEntries(Object.keys(EXPECTED.functional).map((name) => [name, 0]));
  const removedSeqids = new Set();

  await withAtomicGzip(destination, async (writer) => {
    for await (const line of linesFromGzip(source)) {
      if (header === null) {
        header = line.split('\t');
        indexes = Object.fromEntries(header.map((name, index) => [name, index]));
        for (const required of ['protein_id', 'release_seqid_original', ...Object.keys(EXPECTED.functional)]) {
          assert(indexes[required] !== undefined, `functional table lacks required column ${required}`);
        }
        writer.write(`${line}\n`);
        continue;
      }
      if (line === '') continue;
      sourceRows += 1;
      const fields = line.split('\t');
      assert(fields.length === header.length,
        `functional row ${sourceRows} has ${fields.length} fields; expected ${header.length}`);
      const proteinId = fields[indexes.protein_id];
      const originalSeqid = fields[indexes.release_seqid_original];
      const removedByContig = deletedContigs.has(originalSeqid);
      const retainedByProtein = retainedRepresentativeIds.has(proteinId);
      assert(removedByContig === !retainedByProtein,
        `functional row filter disagreement for ${proteinId} on ${originalSeqid}`);
      if (!retainedByProtein) {
        removedRows += 1;
        removedSeqids.add(originalSeqid);
        continue;
      }
      assert(!retainedProteinIds.has(proteinId), `duplicate functional protein_id ${proteinId}`);
      retainedProteinIds.add(proteinId);
      retainedRows += 1;
      for (const metric of Object.keys(counts)) {
        if (nonempty(fields[indexes[metric]])) counts[metric] += 1;
      }
      writer.write(`${line}\n`);
    }
  });

  assert(retainedProteinIds.size === retainedRepresentativeIds.size,
    `functional retained proteins ${retainedProteinIds.size} != representative proteins ${retainedRepresentativeIds.size}`);
  assert(removedSeqids.size === deletedContigs.size,
    `functional rows span ${removedSeqids.size}/${deletedContigs.size} deleted contigs`);
  return { sourceRows, retainedRows, removedRows, retainedProteinIds, counts };
}

async function main() {
  await mkdir(OUT_DIR, { recursive: true });

  const deletionText = await readFile(SOURCES.deletionList, 'utf8');
  const deletedContigs = new Set(
    deletionText.split(/\r?\n/).map((line) => line.trim()).filter((line) => line && !line.startsWith('#')),
  );
  assert(deletedContigs.size === EXPECTED.deletedContigs,
    `deletion list has ${deletedContigs.size} unique contigs; expected ${EXPECTED.deletedContigs}`);

  const genome = await filterFasta(
    SOURCES.genome,
    OUTPUTS.genome,
    (id) => !deletedContigs.has(id),
    true,
  );
  genome.lengthById = new Map();
  {
    // The lengths array follows retained FASTA order; recover an ID->length map in one lightweight pass.
    let currentId = null;
    let currentLength = 0;
    for await (const line of linesFromGzip(OUTPUTS.genome)) {
      if (line.startsWith('>')) {
        if (currentId !== null) genome.lengthById.set(currentId, currentLength);
        currentId = fastaId(line);
        currentLength = 0;
      } else {
        currentLength += line.trim().length;
      }
    }
    if (currentId !== null) genome.lengthById.set(currentId, currentLength);
  }
  const assemblyBp = genome.lengths.reduce((sum, n) => sum + n, 0);
  const n50 = nxStats(genome.lengths, 0.5);
  const n90 = nxStats(genome.lengths, 0.9);
  assert(genome.removedIds.size === EXPECTED.deletedContigs,
    `genome removed ${genome.removedIds.size} contigs; expected ${EXPECTED.deletedContigs}`);
  for (const id of deletedContigs) assert(genome.removedIds.has(id), `deletion-list contig absent from source genome: ${id}`);
  assert(genome.retainedRecords === EXPECTED.contigs,
    `retained contigs ${genome.retainedRecords} != ${EXPECTED.contigs}`);
  assert(assemblyBp === EXPECTED.assemblyBp, `assembly bp ${assemblyBp} != ${EXPECTED.assemblyBp}`);
  assert(n50.nx === EXPECTED.n50 && n50.lx === EXPECTED.l50,
    `N50/L50 ${n50.nx}/${n50.lx} != ${EXPECTED.n50}/${EXPECTED.l50}`);
  assert(n90.nx === EXPECTED.n90 && n90.lx === EXPECTED.l90,
    `N90/L90 ${n90.nx}/${n90.lx} != ${EXPECTED.n90}/${EXPECTED.l90}`);

  const gff = await filterGff(SOURCES.gff, OUTPUTS.gff, deletedContigs);
  assert(gff.genes.size === EXPECTED.genes, `retained genes ${gff.genes.size} != ${EXPECTED.genes}`);
  assert(gff.transcripts.size === EXPECTED.transcripts,
    `retained transcripts ${gff.transcripts.size} != ${EXPECTED.transcripts}`);
  assert(gff.keptFeatureCounts.gene === EXPECTED.genes, 'GFF gene feature count differs from unique gene IDs');
  assert(gff.keptFeatureCounts.mRNA === EXPECTED.transcripts, 'GFF mRNA feature count differs from unique transcript IDs');

  const cdsAll = await filterFasta(SOURCES.cdsAll, OUTPUTS.cdsAll, (id) => gff.transcripts.has(id));
  const proteinAll = await filterFasta(SOURCES.proteinAll, OUTPUTS.proteinAll, (id) => gff.transcripts.has(id));
  assert(cdsAll.retainedRecords === EXPECTED.transcripts,
    `all-isoform CDS count ${cdsAll.retainedRecords} != ${EXPECTED.transcripts}`);
  assert(proteinAll.retainedRecords === EXPECTED.transcripts,
    `all-isoform protein count ${proteinAll.retainedRecords} != ${EXPECTED.transcripts}`);
  assert(cdsAll.retainedIds.size === proteinAll.retainedIds.size &&
    [...cdsAll.retainedIds].every((id) => proteinAll.retainedIds.has(id)),
  'all-isoform CDS and protein ID sets differ');
  assert(cdsAll.retainedIds.size === gff.transcripts.size &&
    [...gff.transcripts].every((id) => cdsAll.retainedIds.has(id)),
  'GFF transcript IDs differ from all-isoform CDS/protein IDs');

  const cdsRep = await filterFasta(SOURCES.cdsRep, OUTPUTS.cdsRep, (id) => gff.transcripts.has(id));
  const proteinRep = await filterFasta(SOURCES.proteinRep, OUTPUTS.proteinRep, (id) => gff.transcripts.has(id));
  assert(cdsRep.retainedRecords === EXPECTED.representatives,
    `representative CDS count ${cdsRep.retainedRecords} != ${EXPECTED.representatives}`);
  assert(proteinRep.retainedRecords === EXPECTED.representatives,
    `representative protein count ${proteinRep.retainedRecords} != ${EXPECTED.representatives}`);
  assert(cdsRep.retainedIds.size === proteinRep.retainedIds.size &&
    [...cdsRep.retainedIds].every((id) => proteinRep.retainedIds.has(id)),
  'representative CDS and protein ID sets differ');
  const representedGenes = new Set([...proteinRep.retainedIds].map((id) => gff.transcriptToGene.get(id)));
  assert(!representedGenes.has(undefined), 'representative protein lacks GFF transcript-to-gene mapping');
  assert(representedGenes.size === EXPECTED.genes,
    `representative sequences cover ${representedGenes.size} genes; expected ${EXPECTED.genes}`);

  const mapping = await filterMapping(SOURCES.mapping, OUTPUTS.mapping, deletedContigs, genome);
  assert(mapping.retainedRows === EXPECTED.contigs,
    `retained mapping rows ${mapping.retainedRows} != ${EXPECTED.contigs}`);

  const functional = await filterFunctional(
    SOURCES.functional,
    OUTPUTS.functional,
    proteinRep.retainedIds,
    deletedContigs,
  );
  assert(functional.retainedRows === EXPECTED.genes,
    `retained functional rows ${functional.retainedRows} != ${EXPECTED.genes}`);
  for (const [metric, expected] of Object.entries(EXPECTED.functional)) {
    assert(functional.counts[metric] === expected,
      `functional ${metric} count ${functional.counts[metric]} != ${expected}`);
  }

  const sourceHashes = {};
  for (const [name, path] of Object.entries(SOURCES)) {
    sourceHashes[name] = {
      path,
      bytes: (await stat(path)).size,
      sha256: await sha256File(path),
    };
  }
  const outputHashes = {};
  for (const [name, path] of Object.entries(OUTPUTS)) {
    outputHashes[name] = {
      file: path.slice(OUT_DIR.length + 1),
      bytes: (await stat(path)).size,
      sha256: await sha256File(path),
    };
  }

  const provenance = {
    schema_version: '1.0',
    formal_sample_id: 'NJ151201',
    species: 'Ephydatia fluviatilis',
    legacy_project_alias: 'S04',
    decision: {
      date: '2026-10-13',
      authority: 'author-confirmed final release decision recorded on 2026-10-13',
      statement_zh: '删除21条支持为非宿主的候选contig作为定稿。',
      action: 'remove every contig listed in E_fluviatilis_removed_21_contigs.txt and all dependent annotation records',
    },
    invariants: {
      source_files_modified: false,
      internal_contig_gene_transcript_ids_renamed: false,
      gzip_level: 9,
      gzip_mtime: 0,
      line_endings: 'LF',
    },
    filters: {
      genome: 'drop FASTA records whose first header token is in the 21-contig deletion set',
      gff3: 'drop all feature rows whose seqid is in the deletion set; retain directives/comments',
      all_isoform_cds_and_protein: 'retain FASTA records whose IDs equal retained GFF3 mRNA IDs',
      representative_cds_and_protein: 'retain source representative FASTA records whose IDs equal retained GFF3 mRNA IDs',
      contig_mapping: 'drop rows whose original_id is in the deletion set',
      functional_table: 'retain rows whose protein_id is in the retained representative protein ID set; verify release_seqid_original is not deleted',
    },
    deleted_contigs: [...deletedContigs],
    validation: {
      assembly: {
        contigs: genome.retainedRecords,
        total_bp: assemblyBp,
        N50_bp: n50.nx,
        L50: n50.lx,
        N90_bp: n90.nx,
        L90: n90.lx,
      },
      annotation: {
        genes: gff.genes.size,
        mRNA: gff.transcripts.size,
        all_isoform_CDS: cdsAll.retainedRecords,
        all_isoform_proteins: proteinAll.retainedRecords,
        representative_CDS: cdsRep.retainedRecords,
        representative_proteins: proteinRep.retainedRecords,
        all_isoform_CDS_protein_ID_sets_equal: true,
        GFF_mRNA_CDS_protein_ID_sets_equal: true,
        representative_CDS_protein_ID_sets_equal: true,
        exactly_one_representative_per_gene: true,
      },
      functional_counts: functional.counts,
      deletion_coverage: {
        deletion_list_unique_contigs: deletedContigs.size,
        genome_removed_contigs: genome.removedIds.size,
        GFF_removed_genes: gff.removedGenes.size,
        GFF_removed_mRNA: gff.removedTranscripts.size,
        mapping_removed_rows: mapping.removedRows,
        functional_removed_rows: functional.removedRows,
      },
      status: 'PASS',
    },
    caveats: [
      'BUSCO, Merqury, RepeatMasker, repeat landscape, mapping/coverage, and contamination-screen metrics were not recomputed here.',
      'Pre-deletion whole-assembly validation metrics must not be represented as measurements of this 21-contig-deleted assembly.',
      'Internal legacy IDs are intentionally retained to preserve exact cross-file links; formal sample ID NJ151201 is used in release filenames and provenance.',
    ],
    sources: sourceHashes,
    generated_outputs: outputHashes,
  };

  const provenancePath = join(OUT_DIR, 'reconstruction_provenance.json');
  await writeTextAtomic(provenancePath, `${JSON.stringify(provenance, null, 2)}\n`);

  const functionalRows = Object.entries(functional.counts)
    .map(([metric, value]) => `| ${metric} | ${value.toLocaleString('en-US')} |`)
    .join('\n');
  const readme = `# NJ151201 final reconstruction\n\n` +
    `Species: *Ephydatia fluviatilis*  \nFormal sample ID: \`NJ151201\`  \nLegacy project alias: \`S04\`\n\n` +
    `This directory is a deterministic derivative built after the author confirmed on 2026-10-13 that deleting the 21 supported non-host candidate contigs is the final decision. The source files were read only and were not modified. Output filenames use the formal sample ID; existing internal contig, gene, and transcript identifiers are retained so that cross-file references remain exact.\n\n` +
    `## Applied rule\n\n` +
    `The genome drops the 21 FASTA records named in \`E_fluviatilis_removed_21_contigs.txt\`. GFF3 rows on those seqids are removed. All-isoform and representative CDS/protein FASTA files are then filtered by the retained GFF3 mRNA ID set. The contig mapping is filtered by \`original_id\`; the functional/model table is filtered by the retained representative protein ID set and independently checked against \`release_seqid_original\`.\n\n` +
    `## Verified final counts\n\n` +
    `- Assembly: 2,097 contigs; 267,354,623 bp; N50 237,588 bp; L50 315; N90 62,222 bp; L90 1,151.\n` +
    `- Structural annotation: 29,180 genes; 33,684 mRNAs; 33,684 all-isoform CDS/protein records.\n` +
    `- Representative set: 29,180 CDS and 29,180 proteins, exactly one retained representative per gene.\n` +
    `- ID validation: retained GFF3 mRNA, all-isoform CDS, and all-isoform protein ID sets are identical; representative CDS/protein ID sets are identical.\n\n` +
    `| Functional metric | Retained genes |\n|---|---:|\n${functionalRows}\n\n` +
    `## Files\n\n` +
    `- \`NJ151201.genome.fna.gz\`: filtered final genome.\n` +
    `- \`NJ151201.annotation.original_contig_ids.gff3.gz\`: filtered standardized GFF3 using original contig IDs.\n` +
    `- \`NJ151201.CDS.all_isoforms.fna.gz\` and \`NJ151201.proteins.all_isoforms.faa.gz\`: retained all-isoform sequences.\n` +
    `- \`NJ151201.CDS.representative.fna.gz\` and \`NJ151201.proteins.representative.faa.gz\`: retained representative sequences.\n` +
    `- \`NJ151201.contig_id_mapping.tsv.gz\`: filtered legacy short-ID/original-ID map.\n` +
    `- \`NJ151201.functional_and_model_annotations.tsv.gz\`: filtered representative-gene functional/model table.\n` +
    `- \`reconstruction_provenance.json\`: source hashes, rule, deletion set, validations, and output hashes.\n` +
    `- \`SHA256SUMS.txt\`: checksums for every deliverable except the checksum file itself.\n\n` +
    `## Important validation boundary\n\n` +
    `No BUSCO, Merqury, RepeatMasker/repeat-landscape, read-mapping/coverage, or contamination-screen result was regenerated. Earlier whole-assembly values describe the pre-deletion assembly and must be recomputed before being attributed to this final derivative. No placeholder or inferred QC result is included.\n`;
  const readmePath = join(OUT_DIR, 'README.md');
  await writeTextAtomic(readmePath, readme);

  const checksumTargets = [
    ...Object.values(OUTPUTS),
    provenancePath,
    readmePath,
  ].sort((a, b) => a.localeCompare(b));
  const checksumLines = [];
  for (const path of checksumTargets) {
    checksumLines.push(`${await sha256File(path)}  ${path.slice(OUT_DIR.length + 1)}`);
  }
  await writeTextAtomic(join(OUT_DIR, 'SHA256SUMS.txt'), `${checksumLines.join('\n')}\n`);

  console.log(JSON.stringify({
    output_directory: OUT_DIR,
    validation_status: 'PASS',
    assembly: provenance.validation.assembly,
    annotation: provenance.validation.annotation,
    functional_counts: functional.counts,
    deletion_coverage: provenance.validation.deletion_coverage,
    generated_outputs: outputHashes,
  }, null, 2));
}

await main();
