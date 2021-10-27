#!/usr/bin/env python

import rulegen.xyz
import rulegen.zstruct
import rulegen.stats
import rulegen.convert

import argparse
import os
import sys


def output_gml_files():
    outdir = "scratch/gmls"
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    _, _, filenames = next(os.walk("scratch/stringfiles"))
    print(filenames)
    for fn in filenames:
        if not fn.startswith("string"):
            continue

        strid = fn[-4:]
        reaction_path = rulegen.xyz.read_stringfile(f"scratch/stringfiles/{fn}")
        gml_str = rulegen.xyz.reaction_to_gml(reaction_path.nodes[0], reaction_path.nodes[-1])
        print(f"{outdir}/{strid}.gml")
        with open(f"{outdir}/{strid}.gml", "w") as f:
            f.write(gml_str)


class RuleGen:
    def __init__(self):
        parser = argparse.ArgumentParser(
            description='Generates zstruct rules',
            usage='''rulegen <command> [<args>]

            The possible commands are
               collect    Runs zstruct and collects the stringfiles
               convert    Converts collected stringfiles into desired format
               stat       Prints relevant stats on collected stringfiles
            ''')
        parser.add_argument('command', help='Subcommand to run')
        # parse_args defaults to [1:] for args, but you need to
        # exclude the rest of the args too, or validation will fail
        args = parser.parse_args(sys.argv[1:2])
        if not hasattr(self, args.command):
            print('Unrecognized command')
            parser.print_help()
            exit(1)
        # use dispatch pattern to invoke method with same name
        getattr(self, args.command)()

    @staticmethod
    def collect():
        parser = argparse.ArgumentParser(
            description='Runs zstruct and collects the stringfiles')
        parser.add_argument('input', help="input file to parse")
        parser.add_argument('--max-count', type=int, help="max number of isomers to apply ssm on")
        args = parser.parse_args(sys.argv[2:])
        isomer_count = rulegen.zstruct.generate_isomers(args.input)
        if args.max_count:
            isomer_count = args.max_count
        rulegen.zstruct.run_ssm(isomer_count)
        output_gml_files()

    @staticmethod
    def stat():
        parser = argparse.ArgumentParser(
            description='Print relevant stats on collected stringfiles')
        parser.add_argument('input', help="stringfiles folder")
        args = parser.parse_args(sys.argv[2:])
        rulegen.stats.aggregate_stringfiles(args.input)

    @staticmethod
    def convert():
        parser = argparse.ArgumentParser(
            description='Converts stringfiles into a chosen format')
        parser.add_argument('input', help="stringfiles folder")
        parser.add_argument('-f', "--format", default="rxn", help="[gml|rxn]")
        parser.add_argument('-o', "--output", default=None, help="output directory")
        args = parser.parse_args(sys.argv[2:])
        if args.format == "rxn":
            rulegen.convert.stringfiles2rxn(args.input, args.output)
        elif args.format == "gml":
            output_gml_files()
        else:
            raise Exception("Unknown format.")


def main():
    RuleGen()


if __name__ == "__main__":
    main()
