#!/usr/bin/env python3
"""
Retorna: nome do processador, núcleos (lógicos e físicos) e memória RAM total.
Requer: py-cpuinfo, psutil
"""

from typing import Optional, Dict, Any
import json

def format_bytes(n: int) -> str:
    """Converte bytes para string legível (GB, MB, etc.)."""
    step = 1024.0
    units = ["B", "KB", "MB", "GB", "TB", "PB"]
    i = 0
    value = float(n)
    while value >= step and i < len(units) - 1:
        value /= step
        i += 1
    return f"{value:.2f} {units[i]}"

def get_cpu_name() -> str:
    """Tenta obter o nome do processador usando py-cpuinfo; faz fallback com platform."""
    try:
        import cpuinfo
        info = cpuinfo.get_cpu_info()
        # chaves comuns retornadas por py-cpuinfo
        for key in ("brand_raw", "brand", "hz_advertised_friendly", "hz_advertised"):
            val = info.get(key)
            if val:
                return str(val).strip()
    except Exception:
        pass

    # fallback simples
    import platform
    cp = platform.processor() or ""
    if cp:
        return cp
    # última tentativa usando uname
    uname = platform.uname()
    for attr in (uname.processor, uname.machine, uname.system):
        if attr:
            return str(attr)
    return "Desconhecido"

def get_core_counts() -> Dict[str, Optional[int]]:
    """
    Retorna dicionário com:
      - logical: número de núcleos lógicos (total)
      - physical: número de núcleos físicos (se disponível)
    """
    import psutil
    import os
    
    logical = psutil.cpu_count(logical=True) or os.cpu_count()
    physical = psutil.cpu_count(logical=False)  # pode retornar None em alguns ambientes
    return {"logical": logical, "physical": physical}

def get_total_ram() -> int:
    """Retorna RAM total em bytes usando psutil."""
    import psutil
    vm = psutil.virtual_memory()
    return int(vm.total)

def system_summary() -> Dict[str, Any]:
    name = get_cpu_name()
    cores = get_core_counts()
    ram_bytes = get_total_ram()

    return {
        "cpu_name": name,
        "cores_logical": cores["logical"],
        "cores_physical": cores["physical"],
        "ram_total_bytes": ram_bytes
        # "ram_total_human": format_bytes(ram_bytes)
    }

if __name__ == "__main__":
    print(system_summary())
