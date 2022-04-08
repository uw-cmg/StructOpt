"""Contains functionality for creating and using loggers."""

import os
import logging


def initialize_logger(filename, name, level=logging.INFO, disable_output=False):
    """Initalizes a logger"""
    logger = logging.getLogger(name)
    logger.setLevel(level)
    if not disable_output:
        if os.path.exists(filename):
            raise RuntimeError("logger filename '{}' already exists".format(filename))
        handler = logging.FileHandler(filename)
        formatter = logging.Formatter("%(asctime)s : %(levelname)s : %(message)s")
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    else:
        return_none = lambda *args, **kwargs: None
        logger.debug = return_none
        logger.info = return_none
        logger.warning = return_none
        logger.error = return_none
        logger.critical = return_none

    return logger


def initialize_logger_for_root(rank, filename="default.log", name="default", level=logging.INFO, disable_output=False):
    if rank == 0:
        logger = initialize_logger(filename=filename, name=name, level=level, disable_output=False)
    else:
        logger = initialize_logger(filename=filename, name=name, level=level, disable_output=True)
    return logger

