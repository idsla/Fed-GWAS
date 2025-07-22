# pipeline/clients/logger_manager.py

import os
import logging
from typing import Optional


class LoggerManager:
    """
    Manages logger instances for federated clients.
    Ensures proper logger initialization and prevents duplicate handlers.
    """
    
    # Class-level cache for logger instances
    _loggers = {}
    
    @classmethod
    def get_logger(cls, client_id: str, log_dir: str) -> logging.Logger:
        """
        Get or create a logger for the specified client.
        
        Args:
            client_id: Unique identifier for the client
            log_dir: Directory where log files should be stored
            
        Returns:
            logging.Logger: Configured logger instance
        """
        logger_name = f"client_{client_id}"
        
        # Check if logger already exists
        if logger_name in cls._loggers:
            return cls._loggers[logger_name]
        
        # Create log directory if it doesn't exist
        os.makedirs(log_dir, exist_ok=True)
        
        # Create new logger
        logger = logging.getLogger(logger_name)
        logger.setLevel(logging.INFO)
        
        # Only add handler if logger doesn't have one
        if not logger.handlers:
            log_path = os.path.join(log_dir, f"client_{client_id}_log.txt")
            file_handler = logging.FileHandler(log_path)
            file_handler.setFormatter(
                logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
            )
            logger.addHandler(file_handler)
        
        # Cache the logger
        cls._loggers[logger_name] = logger
        
        return logger
    
    @classmethod
    def setup_client_logging(cls, client_id: str, log_dir: str, console_logging: bool = False) -> logging.Logger:
        """
        Set up logging for a client with optional console output.
        
        Args:
            client_id: Unique identifier for the client
            log_dir: Directory where log files should be stored
            console_logging: Whether to also log to console
            
        Returns:
            logging.Logger: Configured logger instance
        """
        logger = cls.get_logger(client_id, log_dir)
        
        # Add console handler if requested and not already present
        if console_logging and not any(isinstance(h, logging.StreamHandler) for h in logger.handlers):
            console_handler = logging.StreamHandler()
            console_handler.setFormatter(
                logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
            )
            logger.addHandler(console_handler)
        
        return logger
    
    @classmethod
    def close_logger(cls, client_id: str):
        """
        Close and remove a logger from the cache.
        
        Args:
            client_id: Unique identifier for the client
        """
        logger_name = f"client_{client_id}"
        
        if logger_name in cls._loggers:
            logger = cls._loggers[logger_name]
            
            # Close all handlers
            for handler in logger.handlers[:]:
                handler.close()
                logger.removeHandler(handler)
            
            # Remove from cache
            del cls._loggers[logger_name]