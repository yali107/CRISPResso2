import logging
import sys


def set_logger():
    logging.basicConfig(
        format='%(levelname)-5s @ %(asctime)s:\n\t %(message)s \n',
        datefmt='%a, %d %b %Y %H:%M:%S',
        stream=sys.stderr,
        filemode="w"
    )

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    return logger


logger = set_logger()
