from __future__ import annotations

from Bio import SeqIO
import re
import pandas as pd
import argparse

class InputflError(Exception):
    """Custom exception class for input files.
    Args:
        * Exception ([type]: `str`): Message defining error.
    Returns:
        [type]: `str`: User defined message, or simply error message if 
        no message is defined by the user.
    """

    __module__ = 'builtins'

    def __init__(self, *args):
        if args:
            self.errmessage = args[0]
        else:
            self.errmessage = None

    def __str__(self):
        if self.errmessage:
            return '{0} '.format(self.errmessage)
        else:
            return 'InputflError has been raised'

class _info_parser:
    """Parses all the information required for renaming the fasta file."""

    def __init__(self, csvfl) -> None:
        self.csvfl = csvfl

    @staticmethod
    def _file_parser(fl: str) -> list:
        """Parse through .csv file containing the patterns and their corresponding id and description changes.

        Args:
            * `fl` (str): Input file (.csv format).

        Returns:
            [type]: Three lists;
            * List containing all patterns to find.
            * List containing their corresponding new ids.
            * List containing their corresponding new descriptions.
        """

        df = pd.read_csv(fl)
        dict_df = df.to_dict()
        keys = list(dict_df.keys())

        pattern_key = [keys[0]]
        id_key = [keys[1]]
        desc_key = [keys[2]]

        patterns = [dict_df[y] for y in pattern_key]
        new_ids = [dict_df[y] for y in id_key]
        new_desc = [dict_df[y] for y in desc_key]

        patterns_dict = patterns[0]
        new_ids_dict = new_ids[0]
        new_desc_dict = new_desc[0]

        return patterns_dict, new_ids_dict, new_desc_dict

    def _pattern_parser(self) -> list:
        """Returns a list containing the specified patterns."""
        return list(self._file_parser(fl = self.csvfl)[0].values())

    def _new_ids_parser(self) -> list:
        """Returns a list containing the specified new ids."""
        return list(self._file_parser(fl = self.csvfl)[1].values())

    def _new_desc_parser(self) -> list:
        """Returns a list containing the specified new descriptions."""
        return list(self._file_parser(fl = self.csvfl)[2].values())

    @staticmethod
    def _fasta_info(fl):
        """Calculate and display the number of headers in in input fasta file.

        Args:
            fl ([type]: str): Input fasta file path.
        """

        headers = 0
        for line in fl:
            if line.startswith(">"):
                headers += 1
        print(f"\nThere are {headers} headers in the input fasta.")

class rename_headers(_info_parser):
    def __init__(self, inp :str, out: str, info: bool, change_id: bool, change_desc: bool, csvfl: str) -> None:
        self.inp = inp
        self.out = out
        self.info = info
        self.change_id = change_id
        self.change_desc = change_desc
        super().__init__(csvfl)

    def _record_parser(self, record: str, switch_id: str, switch_description: str, outputfl: str, counter = None) -> None:
        """Header renaming method.

        Args:
            * `record` (str): Element from generator taken from iteration.
            * `switch_id` (str): New fasta header id.
            * `switch_description` (str): New fasta header description.
            * `outputfl` (str): Output fasta with renamed headers.
            * `counter` ([type], optional): Integer var used for enumerating each renamed header.
                                        Useful when many resulting headers have the same name.
                                        Defaults to None.

        Raises:
            ValueError: If counter is specified but value is not of type integer.
        """

        original_header = record.id
        original_description = record.description

        if counter == None:
            counter = ""
        else:
            try:
                int(counter)
            except:
                raise ValueError(f"counter must be of type: int if specified, not of type {type(counter).__name__}.")

        if self.change_id == True:
            record.id = f'{switch_id}{counter}'

        if self.change_desc == True:
            record.description = f'{switch_description}{counter}'

        SeqIO.write(record, outputfl, 'fasta')

        if self.info == True:
            print(f"\nHeader: {original_header} changed to {record.id}. Description changed " 
                ""f"from {original_description} to {record.description}.")

    def _init_parser(self, c: bool, rec: any, iter: int, nid: str, ndesc: str, out: str) -> None:
        """Run _record_parser method by checking boolean condition.

        Args:
            * `c` (bool): Boolean value equal to counter boolean.
            * `rec` (any): Record iterable from generator.
            * `iter` (int): iter number taken from loop iteration in rename()
            * `nid` (str): New id.
            * `ndesc` (str): New description.
            * `out` (str): Output file.
        """

        if c:
            self._record_parser(record = rec, counter = iter, switch_id = nid, 
                    switch_description = ndesc, outputfl = out)

        else:   # no count.
            self._record_parser(record = rec, switch_id = nid, 
                    switch_description = ndesc, outputfl = out)

    def rename(self, count: bool) -> str:
        """Uses regex to find pattern/s on headers.

        Args:
            * `count` (bool): Boolean value of whether or not to write the counter in new names.

        Returns:
            * Path to output fasta file.
        """

        with open(self.inp, "r") as fl, open(self.out, 'w') as outfl:
            tmp = SeqIO.parse(self.inp, 'fasta')   # Use it to list id's.
            records = SeqIO.parse(self.inp, 'fasta')   # Use it for the rest of the operations.

            record_ids = []
            for n in tmp:
                x = n.id
                record_ids.append(x)

            iter = 1
            invalid_patterns = []
            for record in records:  # loop through records.
                for (i, new_id, new_desc) in zip(self._pattern_parser(), self._new_ids_parser(), self._new_desc_parser()):   # nested loop through specified patterns.
                    if not any(i in rec for rec in record_ids):
                        invalid_patterns.append(i)
                        continue

                    if iter > 1: # only when the first pattern in the list is iterated.
                        if re.search(i, record.id):
                            self._init_parser(c = count, rec = record, iter = iter, nid = new_id,
                                            ndesc = new_desc, out = outfl)
                            iter += 1
                        else:
                            continue

                    else:   # Run first pattern seperately.
                        if re.search(i, record.id):
                            self._init_parser(c = count, rec = record, iter = iter, nid = new_id,
                                            ndesc = new_desc, out = outfl)
                            iter += 1
                        else:
                            continue

            # Final checks.
            if iter == 1:
                print("\nNo headers with any of the defined patterns were found.\n")

            else:
                if self.info == True:
                    self._fasta_info(fl)
                    print(f"\n{iter - 1} headers were modified and saved in {self.out}.\n")

            if len(invalid_patterns) > 0:
                str = ', '.join(set(invalid_patterns))
                print(f"The following specified patterns were not detected: {str}\n")

        return outfl

def main(i, cv, o, inf, cnt, id, de):

    if i == None:
        raise InputflError('No fasta file was provided.')
    else:
        inp_fl = i
    if cv == None:
        raise InputflError('No file .csv was provided.')
    else:
        inp_cv = cv

    out_fl = o

    print(inp_fl,inp_cv,out_fl,inf, cnt, id, de)

    rename_headers(inp = inp_fl, out = out_fl, info = inf, change_id = id, change_desc = de, csvfl = inp_cv).rename(count = cnt)


def parse_args(msg):
    parser = argparse.ArgumentParser(description = msg, formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", help = "Input fasta file.")
    parser.add_argument("-cv", help = "Input csv file containing the patterns to look for" 
                        "with their corresponding id's and descriptions to use for renaming.")
    parser.add_argument("-o", default = "output.fasta", type = str, 
                        help = "Optional argument: Name of output fasta file. Default is output.fasta")
    parser.add_argument("-inf", default = False,
                        help = "Optional boolean argument: Display info of header renaming. Default is False.")
    parser.add_argument("-cnt", default = False,
                        help = "Optional boolean argument: Add a counter when naming the new headers. Default option is False.")
    parser.add_argument("-id", default = True,
                        help = "Optional boolean argument: Choose whether to rename header id's. Default is True.")
    parser.add_argument("-de", default = True,
                        help = "Optional boolean argument: Choose whether to rename header descriptions. Default is True.")
    return parser.parse_args()

if __name__ == "__main__":
    msg = ("Header renaming for fasta files.\n\nThis program allows you to rename fasta file headers "
    "using regular expression. To use, please specify the input file path with the -i option,\nthe name of the csv containing"
    " all the required information for renaming and the name of the output file"
    " with the -o option.\nIf no name is specified, the default name will be output.fasta\n"
    "\nThe .csv file should have the following format:\n\n"
    "header_pattern\tnew_id\tnew_description\n"
    "pattern1\tnew_id1\tnew_description1\n"
    "pattern2\tnew_id2\tnew_description2\n\n"
    "\t\t. . .")

    args = parse_args(msg = msg)
    arguments = vars(args)

    i = arguments.get('i')
    cv = arguments.get('cv')
    o = arguments.get('o')
    
    def bool_parse(var):
    
        if type(var) == bool:
            return var
        else:
            if var in ["true", "True", "1"]:
                return True
            elif var in ["false", "False", "0"]:
                return False
            else:
                raise TypeError(f"{var} must be true, True, 1, False, false or 0.")

    inf = bool_parse(arguments.get('inf'))
    cnt = bool_parse(arguments.get('cnt'))
    id = bool_parse(arguments.get('id'))
    de = bool_parse(arguments.get('de'))
    main(i = i, cv = cv, o = o, inf = inf, cnt = cnt, id = id, de = de)