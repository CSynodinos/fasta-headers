#!/usr/bin/env python3
from __future__ import annotations

from Bio import SeqIO
import re
from collections import defaultdict


class headers:
    def __init__(self, inp :str, out: str, pattern: list, new_id: str, new_desc: str) -> None:
        self.inp = inp
        self.out = out
        self.pattern = pattern
        self.new_id = new_id
        self. new_desc = new_desc

    @staticmethod
    def _logic(record: str, switch_id: str, switch_description: str, outputfl: str, counter = None) -> None:
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

        record.id = f'{switch_id}{counter}'
        record.description = f'{switch_description}{counter}'
        SeqIO.write(record, outputfl, 'fasta')

        print(f"\nHeader: {original_header} changed to {record.id}. Description changed" 
            ""f"from {original_description} to {record.description}.")

    @classmethod
    def _run_logic(cls, c: bool, rec: any, iter: int, nid: str, ndesc: str, out: str) -> None:
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
            cls._logic(record = rec, counter = iter, switch_id = nid, 
                    switch_description = ndesc, outputfl = out)

        else:   # no count.
            cls._logic(record = rec, switch_id = nid, 
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
                for i in self.pattern:   # nested loop through specified patterns.
                    if not any(i in rec for rec in record_ids):
                        invalid_patterns.append(i)
                        continue

                    if iter > 1: # only when the first pattern in the list is iterated.
                        if re.search(i, record.id):
                            self._run_logic(c = count, rec = record, iter = iter, nid = self.new_id,
                                            ndesc = self.new_desc, out = outfl)
                            iter += 1
                        else:
                            continue

                    else:   # Run first pattern seperately.
                        if re.search(i, record.id):
                            self._run_logic(c = count, rec = record, iter = iter, nid = self.new_id,
                                            ndesc = self.new_desc, out = outfl)
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
    patterns = [""]
    headers(inp = inp_fl, out = out_fl, pattern = patterns, new_id = "foo", new_desc = "bar").rename()

if __name__ == "__main__":
    main()