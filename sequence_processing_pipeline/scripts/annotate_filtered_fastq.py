#!/usr/bin/env python
import click
from sequence_processing_pipeline.Commands import annotate_filtered_fastq


@click.command()
@click.argument('original_path', required=True)
@click.argument('stripped_path', required=True)
@click.argument('output_path', required=True)
def main(original_path, stripped_path, output_path):
    """
        Annotates a stripped fastq file based on descriptions from original.
        Additional info appears on stderr.
        Exits w/(1) if an id appears in stripped file but not in original.
    """
    try:
        timings = annotate_filtered_fastq(original_path, stripped_path,
                                          output_path)
    except ValueError as e:
        # explicitly catch known potential errors and ensure they exit w/a
        # return value of 1 after politely stating the cause of the error sans
        # stacktrace.
        print(str(e))
        exit(1)

    print(f"Start time: {timings['start_time']}")
    print(f"Time to load optional descriptions: {timings['time_to_load']}")
    print(f"Total time: {timings['total_time']}")


if __name__ == '__main__':
    main()
