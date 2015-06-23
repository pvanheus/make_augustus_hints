#!/usr/bin/python

# from a GFF3 file, generate an augustus hints file

import argparse
from pvh import gff_utils

AUGUSTUS_TYPES = ('start', 'stop', 'tss', 'tts', 'ass', 'dss', 'exonpart', 'exon',
                  'intronpart', 'intron', 'CDSpart', 'CDS', 'UTRpart', 'UTR', 'irpart',
                  'nonexonpart', 'genicpart')

def gff3_to_hints(in_file, out_file, hint_type='XNT', exons_to_CDS=True, trim_cds=15,
                  minintronlen=41, maxintronlen=350000, priority=4,
                  min_transcript_len=500, skip_predicted_proteins=False):
    group = None
    double_trim = trim_cds * 2 # compute this here to save some extra * operations
    transcripts_to_skip = set()
    for line in in_file:
        if line.startswith('#'):
            continue
        line = line.strip()
        fields = line.split('\t')
        assert len(fields) == 9, "Invalid GFF line: {}\n".format(line)
        (ref, source, seq_type, start, end, score, strand, phase, attr_string) = fields
        start = int(start)
        end = int(end)
        # note: in valid GFF3, start is always < end, see http://www.sequenceontology.org/gff3.shtml
        assert start <= end, "start not <= end in GFF line: {}\n".format(line)
        feature_length = end - start + 1
        attributes = gff_utils.parse_gff_attributes(attr_string)
        if seq_type not in AUGUSTUS_TYPES:
            if seq_type == 'transcript':
                # capture the "group id"
                group = attributes['ID']
                # skip the transcript if it is too short or it is a Refseq predicted protein and 
                # we're skipping those
                if ((skip_predicted_proteins and group.startswith('XP_')) 
                        or ((end - start) < min_transcript_len)):
                    transcripts_to_skip.add(group)
            # skip features that AUGUSTUS doesn't care about
            continue
        elif 'Parent' in attributes and attributes['Parent'] in transcripts_to_skip:
            # skip all features that are children of transcripts we want to skip
            continue
        elif exons_to_CDS and seq_type == 'exon':
            seq_type = 'CDS'
        elif seq_type == 'intron' and (feature_length < minintronlen or feature_length > maxintronlen):
            # skip too short or too long introns
            continue
        if seq_type == 'CDS' and trim_cds:
            # if feature is too short (less than 2 * trim_cds), mark new feature in middle of old feature
            seq_type += 'part'
            if feature_length < double_trim:
                start = end = (start + int(feature_length/2))
            else:
                start += trim_cds
                end -= trim_cds
        if group:
            attributes['grp'] = group
        attributes['priority'] = priority
        source = 'xnt2h'
        attributes['src'] = hint_type
        if 'source' in attributes:
            del attributes['source'] # remove the source= attribute
        out_file.write(gff_utils.gff_string_from_list([ref, source, seq_type, start, end, score, strand, phase, attributes]))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert a GFF3 file into augustus hints format')
    parser.add_argument('--hint_type', '-S', default='XNT', help='Set hint type, as per extrinsic config file')
    parser.add_argument('--min_intron_length', type=int, default=41, help='Introns shorter than this length are discarded')
    parser.add_argument('--max_intron_length', type=int, default=350000, help='Introns longer than this are discarded')
    parser.add_argument('--priority', '-P', type=int, default=4)
    parser.add_argument('--source', default='XNT')
    parser.add_argument('--no_exons_to_CDS', default=False, action='store_true', help="Don't convert exon features to CDS")
    parser.add_argument('--CDSpart_cutoff', '-C', type=int, default=15, help='This many bases are cutoff off the start an end of each CDS to make a CDSpart')
    parser.add_argument('--min_transcript_length', type=int, default=500, help='Minimum total transcript length to accept')
    parser.add_argument('--skip_predicted_proteins', '-N', default=False, action='store_true',
                            help='Skip Refseq proteins whose names start with XP_')
    parser.add_argument('gff3_file', type=argparse.FileType())
    parser.add_argument('hints_file', type=argparse.FileType('w'))

    args = parser.parse_args()
    gff3_to_hints(args.gff3_file, args.hints_file, hint_type=args.hint_type, exons_to_CDS=(not args.no_exons_to_CDS), 
                  trim_cds=args.CDSpart_cutoff, minintronlen=args.min_intron_length, maxintronlen=args.max_intron_length,
                  priority=args.priority, min_transcript_len=args.min_transcript_length, 
                  skip_predicted_proteins=args.skip_predicted_proteins)


