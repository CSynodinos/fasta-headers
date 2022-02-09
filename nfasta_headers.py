#!/usr/bin/env python3
from __future__ import annotations

from Bio import SeqIO
import re
import pandas as pd


class _info_parser:
    """Parses all the information required for renaming the fasta file."""

    def __init__(self, csvfl) -> None:
        self.csvfl = csvfl

    @staticmethod
    def _main_parser(fl: str) -> list:
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
        return list(self._main_parser(fl = self.csvfl)[0].values())

    def _new_ids_parser(self) -> list:
        """Returns a list containing the specified new ids."""
        return list(self._main_parser(fl = self.csvfl)[1].values())

    def _new_desc_parser(self) -> list:
        """Returns a list containing the specified new descriptions."""
        return list(self._main_parser(fl = self.csvfl)[2].values())


class headers(_info_parser):
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
            print(f"\nHeader: {original_header} changed to {record.id}. Description changed" 
                ""f"from {original_description} to {record.description}.")

    @classmethod
    def _exec_parser(cls, c: bool, rec: any, iter: int, nid: str, ndesc: str, out: str) -> None:
        """Run _logic static method by checking boolean condition.

        Args:
            * `c` (bool): Boolean value equal to counter boolean.
            * `rec` (any): Record iterable from generator.
            * `iter` (int): iter number taken from loop iteration in rename()
            * `nid` (str): New id.
            * `ndesc` (str): New description.
            * `out` (str): Output file.
        """

        if c:
            cls._record_parser(record = rec, counter = iter, switch_id = nid, 
                    switch_description = ndesc, outputfl = out)

        else:   # no count.
            cls._record_parser(record = rec, switch_id = nid, 
                    switch_description = ndesc, outputfl = out)

    def rename(self, count: bool) -> str:
        """Uses regex to find pattern/s on headers.

        Args:
            * `count` (bool): Boolean value of whether or not to write the counter in new names.

        Returns:
            * Path to output fasta file.
        """


        with open(self.inp, "r") as fl, open(self.out, 'w') as outfl:
            headers = 0
            for line in fl:
                if line.startswith(">"):
                    headers += 1

            tmp = SeqIO.parse(self.inp, 'fasta')   # Use it to list id's.
            records = SeqIO.parse(self.inp, 'fasta')   # Use it for the rest of the operations.
            print(f"\nThere are {headers} headers in the input fasta.")

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
                            self._exec_parser(c = count, rec = record, iter = iter, nid = new_id,
                                            ndesc = new_desc, out = outfl)
                            iter += 1
                        else:
                            continue

                    else:   # Run first pattern seperately.
                        if re.search(i, record.id):
                            self._exec_parser(c = count, rec = record, iter = iter, nid = new_id,
                                            ndesc = new_desc, out = outfl)
                            iter += 1
                        else:
                            continue

        # Final checks.
        if iter == 1:
            print("\nNo headers with any of the defined patterns were found.\n")

        else:
            print(f"\n{iter - 1} headers were modified and saved in {self.out}.\n")

        if len(invalid_patterns) > 0:
            str = ', '.join(set(invalid_patterns))
            print(f"The following specified patterns were not detected: {str}\n")

        return outfl

def main():
    inp_fl = ""
    out_fl = ""
    headers(inp = inp_fl, out = out_fl, info = True, change_id = True, change_desc = True, csvfl = "patterns.csv").rename(count = True)

if __name__ == "__main__":
    main()