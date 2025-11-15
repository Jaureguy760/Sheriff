"""Structured logging configuration for Sheriff."""

import logging
import sys
import json
from pathlib import Path
from typing import Optional, Literal
from datetime import datetime


class JSONFormatter(logging.Formatter):
    """JSON formatter for structured logging."""

    def format(self, record: logging.LogRecord) -> str:
        """Format log record as JSON."""
        log_data = {
            "timestamp": datetime.fromtimestamp(record.created).isoformat(),
            "level": record.levelname,
            "logger": record.name,
            "message": record.getMessage(),
            "module": record.module,
            "function": record.funcName,
            "line": record.lineno,
        }

        if record.exc_info:
            log_data["exception"] = self.formatException(record.exc_info)

        return json.dumps(log_data)


class ColoredFormatter(logging.Formatter):
    """Colored formatter for console output."""

    COLORS = {
        "DEBUG": "\033[36m",  # Cyan
        "INFO": "\033[32m",  # Green
        "WARNING": "\033[33m",  # Yellow
        "ERROR": "\033[31m",  # Red
        "CRITICAL": "\033[35m",  # Magenta
        "RESET": "\033[0m",
    }

    def format(self, record: logging.LogRecord) -> str:
        """Format log record with colors."""
        color = self.COLORS.get(record.levelname, self.COLORS["RESET"])
        reset = self.COLORS["RESET"]

        # Color the level name
        record.levelname = f"{color}{record.levelname}{reset}"

        return super().format(record)


def setup_logging(
    log_file: Optional[str] = None,
    log_level: Literal["DEBUG", "INFO", "WARNING", "ERROR"] = "INFO",
    log_format: Literal["simple", "detailed", "json"] = "detailed",
    console_level: Optional[str] = None,
) -> logging.Logger:
    """Setup Sheriff logging configuration.

    Args:
        log_file: Optional path to log file
        log_level: Logging level for file output
        log_format: Log format (simple, detailed, json)
        console_level: Console log level (defaults to INFO)

    Returns:
        Configured logger instance
    """
    # Get root sheriff logger
    logger = logging.getLogger("sheriff")
    logger.setLevel(logging.DEBUG)  # Capture everything, filter at handler level

    # Remove existing handlers
    logger.handlers.clear()

    # Console handler (always enabled)
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(console_level or "INFO")

    if log_format == "simple":
        console_formatter = logging.Formatter("%(message)s")
    elif log_format == "detailed":
        console_formatter = ColoredFormatter(
            "%(asctime)s [%(levelname)s] %(name)s: %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
        )
    else:  # json
        console_formatter = JSONFormatter()

    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)

    # File handler (optional)
    if log_file:
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)

        file_handler = logging.FileHandler(log_file, mode="a")
        file_handler.setLevel(getattr(logging, log_level))

        if log_format == "json":
            file_formatter = JSONFormatter()
        else:
            file_formatter = logging.Formatter(
                "%(asctime)s [%(levelname)-8s] %(name)-20s %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
            )

        file_handler.setFormatter(file_formatter)
        logger.addHandler(file_handler)

    # Prevent propagation to root logger
    logger.propagate = False

    return logger


def get_logger(name: str) -> logging.Logger:
    """Get a logger with the given name under sheriff hierarchy.

    Args:
        name: Logger name (e.g., 'cli', 'pipeline', 'filtering')

    Returns:
        Logger instance
    """
    return logging.getLogger(f"sheriff.{name}")


# Convenience function for quick setup
def quick_setup(verbosity: int = 1, log_file: Optional[str] = None) -> logging.Logger:
    """Quick logging setup based on verbosity level.

    Args:
        verbosity: 0=ERROR, 1=INFO, 2=DEBUG
        log_file: Optional log file path

    Returns:
        Configured logger
    """
    level_map = {0: "ERROR", 1: "INFO", 2: "DEBUG"}
    level = level_map.get(verbosity, "INFO")

    return setup_logging(log_file=log_file, log_level=level, console_level=level, log_format="detailed")
