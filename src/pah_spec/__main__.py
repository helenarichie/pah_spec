import argparse
import sys
from ._data import add_data_cli_subcommands

def main(args: argparse.Namespace) -> int:
    fn = args.fn
    return fn(args)

parser = argparse.ArgumentParser(description="execute helpful logic for pah_spec")
subparsers = parser.add_subparsers()
add_data_cli_subcommands(subparsers)

if __name__ == "__main__":
    sys.exit(main(parser.parse_args()))
