import subprocess
from collections import defaultdict


def get_run_info(acc_class):
    esearch = subprocess.Popen(['esearch', '-db', 'sra', '-query', acc_class], stdout=subprocess.PIPE)

    efetch = subprocess.Popen(['efetch', '-format', 'runinfo', '-mode', 'xml'], stdin=esearch.stdout,
                              stdout=subprocess.PIPE)
    esearch.wait()
    xtract = subprocess.check_output(['xtract', '-pattern', 'SraRunInfo', '-element', 'LibraryName', 'CenterName'],
                                     stdin=efetch.stdout)
    print(xtract)


run_info = defaultdict(dict)

print(get_run_info('SRR12352757'))
# for accession in open_file_return_lines('SRR_List_1.txt'):
#     acc = Accession(accession)
#     print(acc.acc)
#     print(get_run_info(acc))
