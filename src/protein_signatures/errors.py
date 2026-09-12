"""Custom exceptions for protein signature analysis."""


class ProteinSignatureError(Exception):
    """Base exception for expected package failures."""


class InputValidationError(ProteinSignatureError):
    """Raised when supplied scientific input fails validation."""


class ConfigurationError(ProteinSignatureError):
    """Raised when a campaign configuration is invalid."""


class ExternalToolError(ProteinSignatureError):
    """Raised when an enabled external tool cannot produce valid evidence."""


class PublicationError(ProteinSignatureError):
    """Raised when an immutable result cannot be published safely."""
