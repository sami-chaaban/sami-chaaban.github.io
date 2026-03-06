from __future__ import annotations

from typing import List, Optional
from pydantic import BaseModel, Field


class ChainsRequest(BaseModel):
    pdbId: Optional[str] = Field(default=None, description="PDB identifier")
    mmcifText: Optional[str] = Field(default=None, description="Inline mmCIF text")


class AnalyzeRequest(BaseModel):
    pdbId: Optional[str] = Field(default=None, description="PDB identifier")
    pdbText: Optional[str] = Field(default=None, description="Inline PDB text")
    mmcifText: Optional[str] = Field(default=None, description="Inline mmCIF text")
    chainA: str = Field(..., description="First chain identifier")
    chainB: str = Field(..., description="Second chain identifier")
    mode: str = Field(default="all", description="Filter mode")
    focusResidue: Optional[str] = Field(
        default=None,
        description="Optional focused residue key (e.g. A:318) for residue-scoped analysis acceleration",
    )


class ExplainRequest(BaseModel):
    report: dict
    images: Optional[List[str]] = None
    notes: Optional[str] = None


class RibbonRequest(BaseModel):
    pdbId: Optional[str] = Field(default=None, description="PDB identifier")
    pdbText: Optional[str] = Field(default=None, description="Inline PDB text")
    mmcifText: Optional[str] = Field(default=None, description="Inline mmCIF text")
    step: Optional[float] = Field(default=0.35, description="Sampling step (Å)")


class ChapiMeshRequest(BaseModel):
    pdbId: Optional[str] = Field(default=None, description="PDB identifier")
    pdbText: Optional[str] = Field(default=None, description="Inline PDB text")
    mmcifText: Optional[str] = Field(default=None, description="Inline mmCIF text")
    representation: str = Field(
        default="bonds",
        description="bonds, bonds-selection, ribbon, or surface",
    )

    mode: str = Field(
        default="COLOUR-BY-CHAIN-AND-DICTIONARY",
        description="Bond mesh mode when representation=bonds",
    )
    againstDarkBackground: bool = Field(default=False, description="Adjust colours for dark backgrounds")
    bondWidth: float = Field(default=0.12, description="Bond width in Angstroms")
    atomRadiusToBondWidthRatio: float = Field(default=1.0, description="Ball-and-stick ratio")
    smoothnessFactor: int = Field(default=2, description="Mesh smoothness factor")
    nonDrawCids: Optional[List[str]] = Field(
        default=None,
        description="Atom selection CIDs to exclude from rendering",
    )
    carbonColor: Optional[str] = Field(
        default=None,
        description="Hex color for bespoke carbon atoms",
    )

    cid: str = Field(default="//", description="Atom selection CID for ribbon/surface")
    colourScheme: str = Field(default="Chain", description="Ribbon/surface color scheme")
    style: str = Field(default="Ribbon", description="Ribbon or MolecularSurface")
    secondaryStructureUsage: int = Field(default=0, description="Secondary structure usage flag")
    splitByChain: bool = Field(default=False, description="Return separate meshes per chain")
    chainIds: Optional[List[str]] = Field(
        default=None,
        description="Optional explicit chain-id allowlist for split-by-chain rendering",
    )
