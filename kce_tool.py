import argparse
import sys

from kcet import CTParserByPhase, PkPkiFilter


class KinaseCancerEmbeddingTool(object):
    """
    This class uses a dispathcer pattern to implement one of a list of subcommands used by
    the KCET (KinaseCancerEmbeddingTool) project.
    """

    def __init__(self):
        parser = argparse.ArgumentParser(
            description='kinase cancer embedding tool',
            usage='''kcet <command> [<args>]

The kcet commands are:
   byphase      clinical trials by phase
   pkilist   get list of all protein kinase inhibitors
   pkpki    extract list of pk pki links
''')
        parser.add_argument('command', help='Subcommand to run')
        # parse_args defaults to [1:] for args, but you need to
        # exclude the rest of the args too, or validation will fail
        if len(sys.argv) < 2:
            parser.print_help()
            exit(1)

        args = parser.parse_args(sys.argv[1:2])
        if not hasattr(self, args.command):
            print('Unrecognized command')
            parser.print_help()
            exit(1)
        # use dispatch pattern to invoke method with same name
        getattr(self, args.command)()

    def byphase(self):
        parser = argparse.ArgumentParser(description='parse and integrate yactp data on kinases')
        parser.add_argument('--amend', action='store_true')
        parser.add_argument('-c', '--clinical_trials', required=True,
                            help='Path to the clinical_trials_by_phase.tsv file.')
        parser.add_argument('-y', '--year', type=int, default=None,
                            help='The start year of a clinical trial a drug (minimum 1993).')
        parser.add_argument('--prefix',
                            type=str, default='KCET',
                            help='The prefix for the outfiles.')
        args = parser.parse_args(sys.argv[2:])
        print('Running byphase, yactp=%s' % args.clinical_trials)
        parser = CTParserByPhase(clinical_trials=args.clinical_trials, year=args.year)
        ## Output all phasaes
        filename = "%s_all_phases_%d.tsv" % (args.prefix, parser.get_year())
        df_allphases = parser.get_all_phases()
        df_allphases.to_csv(filename, sep='\t')
        ## Output phase IV
        filename = "%s_phase_4_%d.tsv" % (args.prefix, parser.get_year())
        df_phase4 = parser.get_phase_4()
        df_phase4.to_csv(filename, sep='\t')


    def pkilist(self):
        parser = argparse.ArgumentParser(
            description='extract list of protein kinase inhibitors')
        parser.add_argument('-o', '--outfilename', default='protein_kinase_inhibitors.txt', help="outfilename")
        args = parser.parse_args(sys.argv[2:])
        dklinks = 'input/drug_kinase_links.tsv'
        pki_set = set()
        with open(dklinks) as f:
            header = next(f)
            if not header.startswith('PKI'):
                raise ValueError("Bad header line of drug_kinase_links.tsv: " + header)
            for line in f:
                fields = line.rstrip().split('\t')
                if len(fields) != 4:
                    raise ValueError("Bad line" + line)
                pki = fields[0]
                pki_set.add(pki)
        pki_list = sorted(list(pki_set))
        fh = open(args.outfilename, 'wt')
        for p in pki_list:
            fh.write(p + "\n")
        fh.close()
        print("[INFO] Wrote %d protein kinase inhibitors to %s" % (len(pki_list), args.outfilename))

    def pkpki(self):
        parser = argparse.ArgumentParser(description='Process PKI/PK data')
        parser.add_argument('--max_multiplicity', type=int, default=5)
        parser.add_argument('--outfilename',  type=str, default='input/drug_kinase_links.tsv')
        args = parser.parse_args(sys.argv[2:])
        print(args.max_multiplicity)
        pkpki = PkPkiFilter()
        pkpki.output_to_file(outfilename=args.outfilename, n_pki_limit=args.max_multiplicity)


if __name__ == '__main__':
    KinaseCancerEmbeddingTool()
