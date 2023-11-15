from __future__ import annotations
from typing import Optional, Any, List
from dataclasses import dataclass, field

@dataclass
class Node:
    data: Any
    children: List[Node] = field(default_factory=list)
    parent: Optional[Node] = None


@dataclass
class BiNode:
    data: Any
    children: List[BiNode] = field(default_factory=list)
    parents: List[BiNode] = field(default_factory=list)