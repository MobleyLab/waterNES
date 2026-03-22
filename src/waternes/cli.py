import argparse
from waternes.pipelines.rbfe import run as rbfe_run

def main():
    parser = argparse.ArgumentParser(prog="waternes")
    subparsers = parser.add_subparsers(dest="method")

    # RBFE
    rbfe_parser = subparsers.add_parser("rbfe")
    rbfe_sub = rbfe_parser.add_subparsers(dest="command")

    submit_parser = rbfe_sub.add_parser("submit")
    submit_parser.add_argument("--config", required=True)

    args = parser.parse_args()

    if args.method == "rbfe" and args.command == "submit":
        rbfe_run.submit(args.config)
    else:
        parser.print_help()
