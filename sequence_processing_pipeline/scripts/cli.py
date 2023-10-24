import click
from sequence_processing_pipeline.Commands import demux as _demux


@click.group()
def cli():
    pass


@cli.command()
@click.option('--id-map', type=click.Path(exists=True), required=True)
@click.option('--infile', type=click.Path(exists=True), required=True)
@click.option('--output', type=click.Path(exists=True), required=True)
@click.option('--encoded-id', type=str, required=True)
@click.option('--threads', type=int, required=True)
def demux(id_map, infile, output, encoded_id, threads):
    _demux(id_map, infile, output, encoded_id, threads)


if __name__ == '__main__':
    cli()
