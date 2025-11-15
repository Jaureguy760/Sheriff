"""Checkpoint format definitions."""

from enum import Enum
from typing import Dict, List, Optional, Any
from datetime import datetime
from pydantic import BaseModel, Field
import hashlib
import json


class PipelineStep(str, Enum):
    """Pipeline step identifiers."""

    VALIDATION = "validation"
    BAM_FILTERING = "bam_filtering"
    KMER_MATCHING = "kmer_matching"
    EDIT_CALLING = "edit_calling"
    UMI_COUNTING = "umi_counting"
    OUTPUT_GENERATION = "output_generation"

    @classmethod
    def get_order(cls) -> List[str]:
        """Return steps in execution order."""
        return [
            cls.VALIDATION,
            cls.BAM_FILTERING,
            cls.KMER_MATCHING,
            cls.EDIT_CALLING,
            cls.UMI_COUNTING,
            cls.OUTPUT_GENERATION,
        ]

    @classmethod
    def next_step(cls, current: str) -> Optional[str]:
        """Get next step after current."""
        order = cls.get_order()
        try:
            idx = order.index(current)
            return order[idx + 1] if idx + 1 < len(order) else None
        except ValueError:
            return None


class CheckpointMetrics(BaseModel):
    """Metrics collected during pipeline execution."""

    reads_processed: int = 0
    reads_filtered: int = 0
    barcoded_reads: int = 0
    edit_sites_found: int = 0
    cells_with_edits: int = 0
    genes_quantified: int = 0
    runtime_seconds: float = 0.0


class CheckpointOutputs(BaseModel):
    """Output files and temporary data."""

    filtered_bam: Optional[str] = None
    edit_sites_file: Optional[str] = None
    umi_counts_file: Optional[str] = None
    gene_counts_file: Optional[str] = None
    temp_files: List[str] = Field(default_factory=list)


class Checkpoint(BaseModel):
    """Sheriff pipeline checkpoint."""

    version: str
    timestamp: datetime = Field(default_factory=datetime.now)
    config_hash: str
    status: str = "in_progress"  # in_progress, completed, failed
    current_step: str
    steps_completed: List[str] = Field(default_factory=list)
    steps_remaining: List[str] = Field(default_factory=list)
    outputs: CheckpointOutputs = Field(default_factory=CheckpointOutputs)
    metrics: CheckpointMetrics = Field(default_factory=CheckpointMetrics)
    error: Optional[str] = None

    def to_json(self) -> str:
        """Serialize to JSON string."""
        data = self.model_dump()
        # Convert datetime to ISO format
        data["timestamp"] = self.timestamp.isoformat()
        return json.dumps(data, indent=2)

    @classmethod
    def from_json(cls, json_str: str) -> "Checkpoint":
        """Deserialize from JSON string."""
        data = json.loads(json_str)
        # Convert ISO format back to datetime
        if isinstance(data["timestamp"], str):
            data["timestamp"] = datetime.fromisoformat(data["timestamp"])
        return cls(**data)

    def get_progress_percent(self) -> int:
        """Calculate completion percentage."""
        total_steps = len(PipelineStep.get_order())
        completed = len(self.steps_completed)
        return int((completed / total_steps) * 100)

    @staticmethod
    def compute_config_hash(config: Dict[str, Any]) -> str:
        """Compute hash of configuration for compatibility checking."""
        # Create deterministic string representation
        config_str = json.dumps(config, sort_keys=True)
        return hashlib.sha256(config_str.encode()).hexdigest()[:16]
