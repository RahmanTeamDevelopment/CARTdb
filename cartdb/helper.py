from __future__ import division

from operator import itemgetter
import pysam
import datetime
import os

allowed_chroms = map(str, range(1, 24)) + ['X', 'Y', 'MT']

class TranscriptDBWriter(object):
    """Class for creating new transcript database"""

    def __init__(self, fn, source='', build='', columns=[]):
        """Constructor of the TranscriptDBWriter class"""

        self._fn = fn
        self._source = source
        self._build = build
        self._columns = [x.lower() for x in columns]
        self._records = {c: [] for c in allowed_chroms}
        self.idx_chrom = self._columns.index('chrom')
        self.idx_start = self._columns.index('start')
        self.idx_end = self._columns.index('end')

    def add(self, transcript):
        """Add transcript to DB"""

        record = []
        for c in self._columns:
            if c in ['exons', 'cdna_exons']:
                record.append(','.join([str(e.start) + '-' + str(e.end) for e in getattr(transcript, c.lower())]))
            elif c in ['start', 'end', 'coding_start', 'coding_end', 'cdna_coding_start', 'cdna_coding_end']:
                record.append(int(getattr(transcript, c.lower())))
            else:
                record.append(str(getattr(transcript, c.lower())))
        self._records[transcript.chrom].append(record)

    def _sort_records(self):
        """Sort records by chrom, start, end"""

        idx_start = self._columns.index('start')
        idx_end = self._columns.index('end')
        for c in allowed_chroms:
            if c in self._records:
                self._records[c] = sorted(self._records[c], key=itemgetter(idx_start, idx_end))

    def _index_with_tabix(self):
        """Compress and index output file by Tabix"""

        pysam.tabix_compress(self._fn + '_tmp', self._fn + '.gz', force=True)
        pysam.tabix_index(self._fn + '.gz', seq_col=4, start_col=6, end_col=7, meta_char='#', force=True)


    def finalize(self, options):
        """Write to file, compress and index, clean up"""

        # Sort records by CHROM, START, END
        self._sort_records()

        # Initialize file and write header
        out = open(self._fn + '_tmp', 'w')
        out.write('#createdby: ' + self._source + '\n')
        out.write('#date: ' + str(datetime.datetime.today()).split()[0] + '\n')
        out.write('#build: ' + self._build + '\n')
        out.write('#ncbi_source_db: ' + options.ncbi + '\n')
        out.write('#ucsc_source_db: ' + options.ucsc + '\n')
        out.write('#hgnc_biomart_file: ' + options.hgnc + '\n')

        # Write records to file
        for c in allowed_chroms:
            if c in self._records:
                for record in self._records[c]:
                    record = map(str, record)

                    old_record = [
                        record[0],
                        record[2],
                        record[1],
                        record[3],
                        record[5],
                        '1' if record[4] == '+' else '-1',
                        record[6],
                        record[7],
                        record[-2],
                        str(int(record[9]) + 1),
                        str(int(record[10]) + 1)
                    ]

                    for e in record[8].split(','):
                        [start, end] = e.split('-')
                        old_record.append(start)
                        old_record.append(end)

                    out.write('\t'.join(old_record) + '\n')

        out.close()

        # Compress and index by Tabix
        self._index_with_tabix()

        # Remove temporary file
        os.remove(self._fn + '_tmp')


def output_gff3(transcript, outfile):
    """Output transcript in GFF3 format"""

    attr = ';'.join(['ID=' + transcript.id, 'HGNCID=' + transcript.hgnc_id, 'GENE_SYMBOL=' + transcript.gene_symbol])
    outfile.write('\t'.join([transcript.chrom, '.', 'transcript', str(transcript.start + 1), str(transcript.end), '.', transcript.strand, '.', attr]) + '\n')

    # Exons
    for i in range(len(transcript.exons)):
        exon = transcript.exons[i]
        exon_id = 'EXON' + transcript.id[4:] + '.' + str(i + 1)
        attr = ';'.join(['ID=' + exon_id, 'Parent=' + transcript.id])
        outfile.write('\t'.join([transcript.chrom, '.', 'exon', str(exon.start + 1), str(exon.end), '.', transcript.strand, '.', attr]) + '\n')

    # CDS
    cds_regs = transcript.cds_regions()
    cdspos = 0
    for i in range(len(cds_regs)):
        cds_reg = cds_regs[i]

        cds_id = 'CDS' + transcript.id[4:] + '.' + str(i + 1)
        attr = ';'.join(['ID=' + cds_id, 'Parent=' + transcript.id])

        if cdspos % 3 == 0:
            phase = 0
        elif cdspos % 3 == 1:
            phase = 2
        else:
            phase = 1

        outfile.write('\t'.join([transcript.chrom, '.', 'CDS', str(cds_reg[0] + 1), str(cds_reg[1]), '.', transcript.strand, str(phase), attr]) + '\n')
        cdspos += cds_reg[1] - cds_reg[0]

