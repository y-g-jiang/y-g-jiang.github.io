var Vs=Object.defineProperty;var Xs=(n,t,e)=>t in n?Vs(n,t,{enumerable:!0,configurable:!0,writable:!0,value:e}):n[t]=e;var bt=(n,t,e)=>Xs(n,typeof t!="symbol"?t+"":t,e);function zs(n){const[[t,e,i],[r,s,o],[a,l,u]]=n,h=t*(s*u-o*l)-e*(r*u-o*a)+i*(r*l-s*a);if(Math.abs(h)<1e-10)return null;const c=1/h;return[[(s*u-o*l)*c,(i*l-e*u)*c,(e*o-i*s)*c],[(o*a-r*u)*c,(t*u-i*a)*c,(i*r-t*o)*c],[(r*l-s*a)*c,(e*a-t*l)*c,(t*s-e*r)*c]]}function Ys(n,t){const[[e,i,r],[s,o,a],[l,u,h]]=n,[c,m,d]=t;return[e*c+i*m+r*d,s*c+o*m+a*d,l*c+u*m+h*d]}class Fn{constructor(t){bt(this,"size");bt(this,"isPowerOfTwo");bt(this,"_real");bt(this,"_imag");bt(this,"_scratch",null);bt(this,"_rev",null);bt(this,"_m",0);bt(this,"_internalFFT",null);bt(this,"_chirpReal",null);bt(this,"_chirpImag",null);bt(this,"_bReal",null);bt(this,"_bImag",null);bt(this,"_hanning",null);bt(this,"_windowSumSq",0);this.size=t,this.isPowerOfTwo=(t&t-1)===0&&t>0,this._real=new Float32Array(t),this._imag=new Float32Array(t),this.isPowerOfTwo?this.initRadix2():this.initBluestein()}initRadix2(){const t=this.size,e=Math.log2(t);this._rev=new Uint32Array(t);for(let i=0;i<t;i++){let r=0,s=i;for(let o=0;o<e;o++)r=r<<1|s&1,s>>>=1;this._rev[i]=r}}initBluestein(){const t=this.size;this._m=Math.pow(2,Math.ceil(Math.log2(2*t-1))),this._internalFFT=new Fn(this._m),this._chirpReal=new Float32Array(t),this._chirpImag=new Float32Array(t);for(let r=0;r<t;r++){const s=-Math.PI*(r*r)/t;this._chirpReal[r]=Math.cos(s),this._chirpImag[r]=Math.sin(s)}const e=new Float32Array(this._m),i=new Float32Array(this._m);for(let r=0;r<t;r++)e[r]=this._chirpReal[r],i[r]=-this._chirpImag[r];for(let r=1;r<t;r++)e[this._m-r]=e[r],i[this._m-r]=i[r];this._internalFFT.transform(e,i),this._bReal=new Float32Array(this._internalFFT._real),this._bImag=new Float32Array(this._internalFFT._imag)}initHanning(){if(this._hanning)return;const t=this.size;this._hanning=new Float32Array(t);let e=0;for(let i=0;i<t;i++){const r=.5*(1-Math.cos(2*Math.PI*i/(t-1)));this._hanning[i]=r,e+=r*r}this._windowSumSq=e}transform(t,e){this.isPowerOfTwo?this.transformRadix2(t,e):this.transformBluestein(t,e)}transformRadix2(t,e){const i=this.size,r=this._rev,s=this._real,o=this._imag;if(t===s)for(let a=0;a<i;a++){const l=r[a];if(a<l){const u=s[a],h=o[a];s[a]=s[l],o[a]=o[l],s[l]=u,o[l]=h}}else for(let a=0;a<i;a++){const l=r[a];s[a]=t[l],o[a]=e?e[l]:0}for(let a=2;a<=i;a*=2){const l=a/2,u=-2*Math.PI/a,h=Math.cos(u),c=Math.sin(u);for(let m=0;m<i;m+=a){let d=1,f=0;for(let p=0;p<l;p++){const g=m+p,y=m+p+l,x=d*s[y]-f*o[y],b=d*o[y]+f*s[y],_=s[g],P=o[g];s[g]=_+x,o[g]=P+b,s[y]=_-x,o[y]=P-b;const S=d*h-f*c,w=d*c+f*h;d=S,f=w}}}}transformBluestein(t,e){const i=this.size,r=this._m,s=this._internalFFT,o=s._real,a=s._imag;o.fill(0),a.fill(0);for(let c=0;c<i;c++){const m=t[c],d=e?e[c]:0,f=this._chirpReal[c],p=this._chirpImag[c];o[c]=m*f-d*p,a[c]=m*p+d*f}s.transformRadix2(o,a);for(let c=0;c<r;c++){const m=s._real[c],d=s._imag[c],f=this._bReal[c],p=this._bImag[c];s._real[c]=m*f-d*p,s._imag[c]=m*p+d*f}const l=s._real,u=s._imag;for(let c=0;c<r;c++)u[c]=-u[c];s.transformRadix2(l,u);const h=1/r;for(let c=0;c<i;c++){const m=s._real[c]*h,d=-s._imag[c]*h,f=this._chirpReal[c],p=this._chirpImag[c];this._real[c]=m*f-d*p,this._imag[c]=m*p+d*f}}calculateSpectrum(t,e,i=!1){const r=this.size;let s=0;for(let h=0;h<r;h++)s+=t[h];const o=s/r;this._scratch||(this._scratch=new Float32Array(r));const a=this._scratch;if(i){this.initHanning();const h=this._hanning;for(let c=0;c<r;c++)a[c]=(t[c]-o)*h[c]}else for(let h=0;h<r;h++)a[h]=t[h]-o;this.transform(a);const l=e.length;let u=1/r;i&&this._windowSumSq>0&&(u=1/this._windowSumSq);for(let h=0;h<l;h++){const c=this._real[h],m=this._imag[h];e[h]+=(c*c+m*m)*u}}calculateSpectrumWindow(t,e,i,r=!1){const s=this.size;let o=0;for(let c=0;c<s;c++)o+=t[e+c];const a=o/s;this._scratch||(this._scratch=new Float32Array(s));const l=this._scratch;if(r){this.initHanning();const c=this._hanning;for(let m=0;m<s;m++)l[m]=(t[e+m]-a)*c[m]}else for(let c=0;c<s;c++)l[c]=t[e+c]-a;this.transform(l);const u=i.length;let h=1/s;r&&this._windowSumSq>0&&(h=1/this._windowSumSq);for(let c=0;c<u;c++){const m=this._real[c],d=this._imag[c];i[c]+=(m*m+d*d)*h}}}const Ws={"Sony ILCE-7RM5":"0.82 -0.2976 -0.0719 -0.4296 1.2053 0.2532 -0.0429 0.1282 0.5774"};let Pi=null;async function Hs(n){return Pi||(Pi=(async()=>{if(typeof window.loadPyodide!="function")throw new Error("Pyodide missing: window.loadPyodide not found.");const t=await window.loadPyodide();return await t.loadPackage("numpy"),t})()),Pi}var js=`#!/usr/bin/env python3
"""Container parser for the Sony ARW6/LLVC3 raw strip.

This stops at packet records and control fields. Coefficient entropy decoding
lives in llvc3_entropy.py; this file is mostly the thing I dump to JSON and diff
against Imaging Edge traces when the packet framing looks suspicious.
"""

from __future__ import annotations

import argparse
import json
import struct
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


RAW_STREAM_OFFSET = 0x200


def u16be(buf: bytes, off: int) -> int:
    return struct.unpack_from(">H", buf, off)[0]


def u16le(buf: bytes, off: int) -> int:
    return struct.unpack_from("<H", buf, off)[0]


def u32le(buf: bytes, off: int) -> int:
    return struct.unpack_from("<I", buf, off)[0]


def u32be(buf: bytes, off: int) -> int:
    return struct.unpack_from(">I", buf, off)[0]


class BitReader:
    """Small MSB-first reader for packet headers."""

    def __init__(self, data: bytes) -> None:
        self.data = data
        self.pos = 0

    def read(self, nbits: int) -> int:
        if nbits < 0:
            raise ValueError("negative bit count")
        out = 0
        for _ in range(nbits):
            if self.pos >= len(self.data) * 8:
                raise EOFError("packet bitstream exhausted")
            byte = self.data[self.pos >> 3]
            bit = (byte >> (7 - (self.pos & 7))) & 1
            out = (out << 1) | bit
            self.pos += 1
        return out


@dataclass
class TiffRawInfo:
    ifd_offset: int
    width: int
    height: int
    bits_per_sample: int
    compression: int
    photometric: int
    strip_offset: int
    strip_byte_count: int
    cfa_pattern: list[int] | None
    black_level_tag_0x7310: list[int] | None
    white_level: int | None
    default_crop_origin: list[int] | None
    default_crop_size: list[int] | None


@dataclass
class LlvcHeader:
    magic: str
    sequence_or_version: int
    coded_width: int
    coded_half_height: int
    logical_height: int
    decoded_bits: int
    component_count: int
    mode: int
    flags_low10: int


@dataclass
class LlvcStreamInfo:
    index: int
    offset: int
    length: int
    header: LlvcHeader
    tile_x: int
    tile_y: int
    tile_width: int
    tile_height: int


@dataclass
class DirectoryEntry:
    group: int
    index: int
    start: int
    length: int
    mode_in_decoder_trace: int | None = None


@dataclass
class PacketRecord:
    index: int
    byte_length: int
    selectors: list[int]
    payload_offset: int


@dataclass
class PacketHeader:
    group: int
    index: int
    stream_offset: int
    directory_length: int
    control_count: int
    extra_count: int
    tag4: int
    reserved4: int
    type2: int
    control_words: int
    block_count: int
    width_marker: int
    reserved8: int
    skipped_u8: list[int]
    header_bits: int
    control_bytes: int
    total_bytes: int
    payload_bytes_from_records: int
    validation: dict[str, bool]
    records: list[PacketRecord]
    first_records: list[PacketRecord]
    last_records: list[PacketRecord]


def parse_tiff_value(data: bytes, bo: str, typ: int, cnt: int, raw: bytes) -> Any:
    fmt_by_type = {1: "B", 3: "H", 4: "I", 8: "h", 9: "i"}
    if typ in fmt_by_type:
        fmt = bo + fmt_by_type[typ] * cnt
        size = struct.calcsize(fmt)
        val = list(struct.unpack(fmt, raw[:size]))
        return val[0] if len(val) == 1 else val
    if typ in (2, 7):
        return list(raw[:cnt])
    return raw.hex(" ")


def iter_tiff_ifds(data: bytes) -> list[tuple[int, dict[int, Any]]]:
    if data[:2] == b"II":
        bo = "<"
    elif data[:2] == b"MM":
        bo = ">"
    else:
        raise ValueError("input is not a TIFF/ARW file")

    def u16(off: int) -> int:
        return struct.unpack_from(bo + "H", data, off)[0]

    def u32(off: int) -> int:
        return struct.unpack_from(bo + "I", data, off)[0]

    type_sizes = {1: 1, 2: 1, 3: 2, 4: 4, 5: 8, 7: 1, 8: 2, 9: 4, 10: 8}
    stack = [u32(4)]
    seen: set[int] = set()
    out: list[tuple[int, dict[int, Any]]] = []
    while stack:
        off = stack.pop()
        if off in seen or off <= 0 or off >= len(data):
            continue
        seen.add(off)
        n = u16(off)
        tags: dict[int, Any] = {}
        for i in range(n):
            ent = off + 2 + i * 12
            tag = u16(ent)
            typ = u16(ent + 2)
            cnt = u32(ent + 4)
            value_area = ent + 8
            size = type_sizes.get(typ, 1) * cnt
            raw = data[value_area : value_area + 4] if size <= 4 else data[u32(value_area) : u32(value_area) + size]
            tags[tag] = parse_tiff_value(data, bo, typ, cnt, raw)
            if tag == 0x014A:
                offsets = tags[tag] if isinstance(tags[tag], list) else [tags[tag]]
                stack.extend(int(x) for x in offsets)
        next_ifd = u32(off + 2 + n * 12)
        if next_ifd:
            stack.append(next_ifd)
        out.append((off, tags))
    return out


def scalar_or_first(value: Any) -> int:
    if isinstance(value, list):
        return int(value[0])
    return int(value)


def find_raw_subifd(path: Path) -> tuple[TiffRawInfo, bytes]:
    data = path.read_bytes()
    for off, tags in iter_tiff_ifds(data):
        if tags.get(0x0103) == 32766 and tags.get(0x0106) == 32803:
            strip_offset = scalar_or_first(tags[0x0111])
            strip_len = scalar_or_first(tags[0x0117])
            info = TiffRawInfo(
                ifd_offset=off,
                width=int(tags[0x0100]),
                height=int(tags[0x0101]),
                bits_per_sample=int(tags[0x0102]),
                compression=int(tags[0x0103]),
                photometric=int(tags[0x0106]),
                strip_offset=strip_offset,
                strip_byte_count=strip_len,
                cfa_pattern=list(tags.get(0x828E, [])) or None,
                black_level_tag_0x7310=list(tags.get(0x7310, [])) if isinstance(tags.get(0x7310), list) else None,
                white_level=int(tags[0xC61D]) if 0xC61D in tags else None,
                default_crop_origin=list(tags.get(0xC61F, [])) if isinstance(tags.get(0xC61F), list) else None,
                default_crop_size=list(tags.get(0xC620, [])) if isinstance(tags.get(0xC620), list) else None,
            )
            return info, data[strip_offset : strip_offset + strip_len]
    raise ValueError("no ARW6 LLVC raw SubIFD found")


def parse_llvc_header(stream: bytes) -> LlvcHeader:
    word_c = u16be(stream, 0x0C)
    word_e = u16be(stream, 0x0E)
    return LlvcHeader(
        magic=stream[:4].decode("ascii", "replace"),
        sequence_or_version=u32le(stream, 0x04),
        coded_width=u16be(stream, 0x08),
        coded_half_height=u16be(stream, 0x0A),
        logical_height=u16be(stream, 0x0A) * 2,
        decoded_bits=(word_c >> 4) & 0x3F,
        component_count=word_e >> 13,
        mode=(word_e >> 10) & 0x03,
        flags_low10=word_e & 0x03FF,
    )


def initial_group_lengths(stream: bytes) -> list[int]:
    return [((u32be(stream, 0x10 + off) >> 4) & 0x0FFFFFF0) for off in (0, 3, 6, 9, 12)]


def llvc_stream_length(stream: bytes) -> int:
    """Return the byte span occupied by one LLVC3 stream."""

    return 0x80 + sum(initial_group_lengths(stream))


def is_plausible_llvc_header(header: LlvcHeader) -> bool:
    return (
        header.magic in {"A000", "0000"}
        and header.coded_width > 0
        and header.coded_half_height > 0
        and header.decoded_bits == 16
        and header.component_count == 3
        and header.mode == 3
    )


def find_llvc_streams(strip: bytes) -> list[LlvcStreamInfo]:
    """Find all LLVC3 streams inside an ARW6 raw strip."""

    streams: list[LlvcStreamInfo] = []
    count = u32le(strip, 0) if len(strip) >= 4 else 0
    if 1 <= count <= 16 and len(strip) >= RAW_STREAM_OFFSET + 0x80:
        for index in range(count):
            entry = 0x08 + index * 0x18
            if entry + 0x18 > len(strip):
                streams = []
                break
            table_offset = u32le(strip, entry)
            tile_x = u32le(strip, entry + 0x08)
            tile_y = u32le(strip, entry + 0x0C)
            tile_width = u32le(strip, entry + 0x10)
            tile_height = u32le(strip, entry + 0x14)
            pos = table_offset if table_offset else RAW_STREAM_OFFSET
            found_pos: int | None = None
            search_end = min(len(strip) - 0x80, pos + 0x1000)
            for candidate in range(pos, search_end + 1, 0x10):
                try:
                    candidate_header = parse_llvc_header(strip[candidate:])
                except Exception:
                    continue
                if is_plausible_llvc_header(candidate_header):
                    found_pos = candidate
                    header = candidate_header
                    break
            if found_pos is None:
                streams = []
                break
            pos = found_pos
            if not is_plausible_llvc_header(header):
                streams = []
                break
            length = llvc_stream_length(strip[pos:])
            if length <= 0x80 or pos + length > len(strip):
                streams = []
                break
            streams.append(
                LlvcStreamInfo(
                    index=index,
                    offset=pos,
                    length=length,
                    header=header,
                    tile_x=tile_x,
                    tile_y=tile_y,
                    tile_width=tile_width or header.coded_width,
                    tile_height=tile_height or header.logical_height,
                )
            )
        if len(streams) == count:
            return streams

    for pos in range(RAW_STREAM_OFFSET, max(RAW_STREAM_OFFSET, len(strip) - 0x80), 0x10):
        try:
            header = parse_llvc_header(strip[pos:])
        except Exception:
            continue
        if not is_plausible_llvc_header(header):
            continue
        try:
            length = llvc_stream_length(strip[pos:])
        except Exception:
            continue
        if length <= 0x80 or pos + length > len(strip):
            continue
        streams.append(
            LlvcStreamInfo(
                index=len(streams),
                offset=pos,
                length=length,
                header=header,
                tile_x=0,
                tile_y=0,
                tile_width=header.coded_width,
                tile_height=header.logical_height,
            )
        )
    return streams


def parse_directory(stream: bytes) -> tuple[int, list[DirectoryEntry], list[dict[str, Any]]]:
    group_lengths = initial_group_lengths(stream)
    pos = 0x30
    base = 0
    entries: list[DirectoryEntry] = []
    groups: list[dict[str, Any]] = []
    for group, group_len in enumerate(group_lengths):
        n_entries = stream[pos] & 0x0F
        values = [(u32be(stream, pos + off) & 0x00FFFFFF) << 4 for off in (0, 3, 6, 9, 12)]
        consumed = 0x10
        if n_entries >= 5:
            values.extend((u32be(stream, pos + 0x10 + off) & 0x00FFFFFF) << 4 for off in (0, 3, 6, 9))
            consumed = 0x20
        local = 0
        for index in range(n_entries):
            entries.append(DirectoryEntry(group=group, index=index, start=0x80 + base + local, length=values[index]))
            local += values[index]
        groups.append(
            {
                "group": group,
                "declared_length": group_len,
                "entry_count": n_entries,
                "entry_lengths": values[:n_entries],
                "entry_sum": local,
                "directory_offset": pos,
                "directory_bytes": consumed,
                "sum_matches_declared": local == group_len,
            }
        )
        base += group_len
        pos += consumed
    return pos, entries, groups


def parse_packet(stream: bytes, entry: DirectoryEntry) -> PacketHeader:
    packet = stream[entry.start : entry.start + entry.length]
    br = BitReader(packet)
    control_count = br.read(16)
    extra_count = br.read(24)
    tag4 = br.read(4)
    reserved4 = br.read(4)
    type2 = br.read(2)
    control_words = br.read(6)
    block_count = br.read(16)
    width_marker = br.read(8)
    reserved8 = br.read(8)
    skipped = [br.read(8) for _ in range(5)]
    control_bytes = (control_count + 1) << 4
    total_bytes = (control_count + 1 + extra_count) << 4
    records: list[PacketRecord] = []
    cursor = control_bytes
    for i in range(block_count):
        byte_len = br.read(16)
        selectors = [br.read(4) for _ in range(type2)]
        records.append(PacketRecord(i, byte_len, selectors, cursor))
        cursor += byte_len
    # Native code consumes the 6-bit field while walking the header, then the
    # validation at 0x1aa717 compares the original 16-bit control_count with:
    # ceil((type2 + 4) * block_count * 4 / 128).
    formula = ((type2 + 4) * block_count * 4 + 0x7F) >> 7
    validation = {
        "directory_length_matches_total": entry.length == total_bytes,
        "tag4_is_4": tag4 == 4,
        "reserved4_is_0": reserved4 == 0,
        "type2_is_1_or_3": type2 in (1, 3),
        "reserved6_is_0": control_words == 0,
        "control_count_formula": control_count == formula,
        "block_count_le_300": block_count <= 300,
        "width_marker_is_0x10": width_marker == 0x10,
        "reserved8_is_0": reserved8 == 0,
        "skipped_bytes_are_0": all(x == 0 for x in skipped),
        "record_payload_within_total": cursor <= total_bytes,
    }
    return PacketHeader(
        group=entry.group,
        index=entry.index,
        stream_offset=entry.start,
        directory_length=entry.length,
        control_count=control_count,
        extra_count=extra_count,
        tag4=tag4,
        reserved4=reserved4,
        type2=type2,
        control_words=control_words,
        block_count=block_count,
        width_marker=width_marker,
        reserved8=reserved8,
        skipped_u8=skipped,
        header_bits=br.pos,
        control_bytes=control_bytes,
        total_bytes=total_bytes,
        payload_bytes_from_records=cursor,
        validation=validation,
        records=records,
        first_records=records[:8],
        last_records=records[-4:],
    )


def derive_metrics(raw_info: TiffRawInfo, strip_len: int, stream_len: int) -> dict[str, Any]:
    samples = raw_info.width * raw_info.height
    return {
        "samples": samples,
        "decoded_u16_bytes": samples * 2,
        "decoded_14bit_packed_bytes": samples * 14 / 8,
        "strip_bytes_per_sample": strip_len / samples,
        "strip_bits_per_sample": strip_len * 8 / samples,
        "stream_bits_per_sample_excluding_0x200_preamble": stream_len * 8 / samples,
        "compression_ratio_vs_u16": (samples * 2) / strip_len,
        "compression_ratio_vs_14bit_packed": (samples * 14 / 8) / strip_len,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("path", nargs="?", default="DSC00089.ARW")
    ap.add_argument("--out", default="out/reverse/llvc3_bitstream_probe.json")
    ns = ap.parse_args()

    raw_info, strip = find_raw_subifd(Path(ns.path))
    stream = strip[RAW_STREAM_OFFSET:]
    header = parse_llvc_header(stream)
    consumed, entries, groups = parse_directory(stream)
    packets = [parse_packet(stream, entry) for entry in entries]
    all_valid = all(all(p.validation.values()) for p in packets)
    type_counts: dict[str, int] = {}
    for p in packets:
        type_counts[str(p.type2)] = type_counts.get(str(p.type2), 0) + 1
    result = {
        "input": str(ns.path),
        "raw_subifd": asdict(raw_info),
        "llvc_header": asdict(header),
        "stream_offset_inside_raw_strip": RAW_STREAM_OFFSET,
        "strip_preamble_first_32": strip[:32].hex(" "),
        "directory": {
            "consumed_bytes_after_stream_header": consumed - 0x10,
            "packet_base_offset": 0x80,
            "groups": groups,
            "entries": [asdict(e) for e in entries],
        },
        "packets": [asdict(p) for p in packets],
        "summary": {
            "packet_count": len(packets),
            "packet_type_counts": type_counts,
            "all_packet_validations_pass": all_valid,
            "total_packet_bytes": sum(p.directory_length for p in packets),
            "payload_bytes_from_records": sum(p.payload_bytes_from_records for p in packets),
            "metrics": derive_metrics(raw_info, len(strip), len(stream)),
        },
    }
    out = Path(ns.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result["summary"], indent=2))


if __name__ == "__main__":
    main()
`,Qs=`#!/usr/bin/env python3
"""Entropy-side notes for the Sony ARW6/LLVC3 stream.

The bit reader and 4-lane coefficient paths here came straight out of Imaging
Edge traces. I kept the code narrow on purpose: first replay a row, then a
packet, then let the higher-level decoder stitch the pieces together.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from llvc3_bitstream_probe import find_llvc_streams, find_raw_subifd


@dataclass
class NativeBitReader:
    """MSB-first 64-bit reader with the same state layout Sony uses.

    The offsets are worth keeping close by:

      +0x00 current big-endian 64-bit word
      +0x08 pointer to the last loaded word
      +0x10 bit position / remaining bits in current word
      +0x14 64-bit words left to load
      +0x1c status

    The pointer is initialized to eight bytes before the record payload. On
    underflow the decoder advances by eight bytes and loads a big-endian word.
    """

    data: bytes
    ptr: int = -8
    cur: int = 0
    bit: int = 0
    words_left: int = 0
    status: int = 0

    @classmethod
    def for_record(cls, packet: bytes, payload_offset: int, byte_length: int) -> "NativeBitReader":
        words = (((byte_length + 7) // 8) + 1) & ~1
        return cls(data=packet, ptr=payload_offset - 8, words_left=words)

    def _load_next_word(self) -> int:
        self.words_left -= 1
        if self.words_left < 0:
            self.status = 2
            return 0
        self.ptr += 8
        chunk = self.data[self.ptr : self.ptr + 8]
        if len(chunk) < 8:
            chunk = chunk + b"\\x00" * (8 - len(chunk))
        self.cur = int.from_bytes(chunk, "big")
        self.bit = 64
        return self.cur

    def read_bits(self, nbits: int) -> int:
        if nbits == 0:
            return 0
        if nbits < 0:
            raise ValueError("negative bit count")
        out = 0
        remaining = nbits
        while remaining > 0:
            if self.bit <= 0:
                self._load_next_word()
                if self.status:
                    return out << remaining
            take = min(remaining, self.bit)
            self.bit -= take
            out = (out << take) | ((self.cur >> self.bit) & ((1 << take) - 1))
            remaining -= take
        return out

    def skip_zero_words_to_one(self, initial_bit_count: int) -> int:
        """Count leading zeros until the next one, in the 0x1a9080 style."""

        total = 0
        while True:
            if self.bit <= 0:
                self._load_next_word()
                if self.status:
                    return total
            if self.cur == 0:
                total += self.bit
                self.bit = 0
                continue
            # Count zeros among currently available MSB-side bits.
            while self.bit > 0:
                self.bit -= 1
                if (self.cur >> self.bit) & 1:
                    return total
                total += 1
                if total >= initial_bit_count and initial_bit_count > 0:
                    return total

    def read_unary_zeros_plus_one(self) -> int:
        zeros = 0
        while not self.status:
            bit = self.read_bits(1)
            if self.status:
                break
            if bit:
                return zeros + 1
            zeros += 1
        return zeros + 1


def split_lanes(packed: int, width: int) -> list[int]:
    mask = (1 << width) - 1
    return [
        (packed >> (3 * width)) & mask,
        (packed >> (2 * width)) & mask,
        (packed >> width) & mask,
        packed & mask,
    ]


def magnitude4(br: NativeBitReader, width: int, shift: int) -> list[int]:
    """Magnitude expansion from 0x1a8b00."""

    if width <= 0:
        raw = [0, 0, 0, 0]
    else:
        raw = split_lanes(br.read_bits(width * 4), width)
    if shift <= 0:
        return raw
    out: list[int] = []
    for x in raw:
        if x <= 0:
            out.append(0)
        else:
            # Visible SIMD shape: ((2*x + 1) << (shift - 1)) - (x & 1).
            out.append(((2 * x + 1) << (shift - 1)) - (x & 1))
    return out


def apply_sign4(br: NativeBitReader, coeffs: Iterable[int]) -> list[int]:
    """Apply one sign bit to each positive lane; see 0x1a8dd0."""

    out: list[int] = []
    for x in coeffs:
        if x > 0:
            bit = br.read_bits(1)
            out.append(x - 2 * x * bit)
        else:
            out.append(x)
    return out


def update_width(br: NativeBitReader, width: int) -> int:
    """Adaptive width update shared by 0x1a9080 and 0x1ac060.

    Prefix structure from traces and branch shape:

      0      keep width
      10 U   increase by unary U, where U is zeros+1 terminated by a one
      11 U   decrease by unary U, clipped at zero
    """

    if br.read_bits(1) == 0 or br.status:
        return width
    if br.read_bits(1) == 0 or br.status:
        return width + br.read_unary_zeros_plus_one()

    # Decrement path: if the run would cross zero, native code consumes only
    # width-1 zeros and leaves the terminating one for the next state; that is
    # the 0x1ac4ff/0x1ac50b saturation branch.
    for zeros in range(max(0, width - 1)):
        if br.read_bits(1):
            return width - (zeros + 1)
    return 0


def read_initial_width(br: NativeBitReader) -> int:
    """Read the row's first adaptive width with the 0x1a9080 zero state."""

    return update_width(br, 0)


def read_zero_run(br: NativeBitReader, remaining: int) -> int:
    """Read the zero-group run used when the adaptive width is zero."""

    if remaining <= 1:
        return remaining
    max_prefix = (remaining - 1).bit_length()
    zeros = 0
    while zeros < max_prefix:
        if br.read_bits(1):
            break
        if br.status:
            return remaining
        zeros += 1
    else:
        return remaining
    base = 1 << zeros
    if base >= remaining:
        return remaining
    extra = br.read_bits(zeros) if zeros else 0
    run = base + extra
    return min(run, remaining)


def packet_row_multiplier(group: int) -> int:
    """Number of output rows produced by one packet record for this scale."""

    if group < 0:
        raise ValueError("negative packet group")
    if group == 4:
        return 8
    if group == 0:
        return 1
    return 1 << (group - 1)


def read_width_after_zero_run(br: NativeBitReader) -> int:
    """Positive width code after a zero run: leading zeros plus the one bit."""

    return br.read_unary_zeros_plus_one()


def alt4(br: NativeBitReader, width: int, shift: int, has_next: bool = True) -> tuple[list[int], int]:
    """Current reconstruction of the alternate 0x1ac060 4-lane path.

    It shares the magnitude transform with 0x1a8b00, then signs only positive
    lanes. The width predictor agrees with the small-width traces I have, but
    wider rows still deserve spot checks.
    """

    coeffs = magnitude4(br, width, shift)
    next_width = width
    if has_next:
        next_width = update_width(br, width)
    return apply_sign4(br, coeffs), next_width


def decode_record_component(
    br: NativeBitReader, groups: int, shift: int = 0
) -> tuple[list[int], list[int], int]:
    """Decode one component from the current bitreader position.

    \`groups\` is the number of 4-lane coefficient groups. The packet selector
    nibble is passed as \`shift\`; the initial adaptive width is read from the
    record payload itself by the 0x1a9080 prefix reader.
    """

    width = read_initial_width(br)
    initial_width = width
    coeffs: list[int] = []
    widths: list[int] = []
    gi = 0
    while gi < groups:
        if br.status or width > 0x13:
            br.status = br.status or 1
            coeffs.extend([0, 0, 0, 0] * (groups - gi))
            widths.extend([width] * (groups - gi))
            break
        if width == 0:
            run = read_zero_run(br, groups - gi)
            coeffs.extend([0, 0, 0, 0] * run)
            widths.extend([0] * run)
            gi += run
            if gi >= groups:
                break
            width = read_width_after_zero_run(br)
            continue
        vals, width = alt4(br, width, shift=shift, has_next=gi + 1 < groups)
        coeffs.extend(vals)
        widths.append(width)
        gi += 1
    return coeffs, widths, initial_width


def load_packet(arw: Path, group: int, index: int, stream_index: int = 0) -> tuple[bytes, dict]:
    from llvc3_bitstream_probe import parse_directory, parse_packet

    raw_info, strip = find_raw_subifd(arw)
    streams = find_llvc_streams(strip)
    if not streams:
        raise ValueError("no LLVC3 stream found in ARW6 raw strip")
    if stream_index < 0 or stream_index >= len(streams):
        raise ValueError(f"stream_index {stream_index} out of range for {len(streams)} LLVC3 streams")
    stream_info = streams[stream_index]
    stream = strip[stream_info.offset : stream_info.offset + stream_info.length]
    _consumed, entries, _groups = parse_directory(stream)
    entry = next(e for e in entries if e.group == group and e.index == index)
    packet_info = parse_packet(stream, entry)
    packet = stream[entry.start : entry.start + entry.length]
    info = json.loads(json.dumps(packet_info, default=lambda o: o.__dict__))
    header = stream_info.header
    info["raw_width"] = raw_info.width
    info["raw_height"] = raw_info.height
    info["work_width"] = header.coded_width
    info["work_height"] = header.logical_height
    info["stream_index"] = stream_index
    info["stream_offset"] = stream_info.offset
    info["stream_length"] = stream_info.length
    info["llvc_header"] = header.__dict__
    return packet, info


def replay_row(arw: Path, group: int, index: int, row: int, groups: int = 16) -> dict[str, object]:
    packet, info = load_packet(arw, group, index)
    rec = info["records"][row]
    coeffs, widths, br, initial_width = decode_type1_row(packet, rec, groups)
    return {
        "packet": {"group": group, "index": index},
        "row": row,
        "record": rec,
        "selector": rec["selectors"][0] if rec["selectors"] else 0,
        "initial_width": initial_width,
        "coeffs": coeffs,
        "next_widths": widths,
        "bitreader": {"ptr": br.ptr, "bit": br.bit, "words_left": br.words_left, "status": br.status},
    }


def decode_type1_row(packet: bytes, rec: dict, groups: int) -> tuple[list[int], list[int], NativeBitReader, int]:
    """Decode one type-1 record into 4-lane signed coefficients."""

    br = NativeBitReader.for_record(packet, rec["payload_offset"], rec["byte_length"])
    shift = rec["selectors"][0] if rec["selectors"] else 0
    coeffs, widths, initial_width = decode_record_component(br, groups, shift)
    return coeffs, widths, br, initial_width


def decode_record_components(
    packet: bytes, rec: dict, groups: int, components: int, row_multiplier: int = 1
) -> tuple[list[list[int]], list[dict[str, int]]]:
    """Decode all components stored in one packet record.

    Type-3 records share one payload bitreader across three component streams.
    The native block decoder calls the same entropy routine three times with
    selector nibbles from the control record, so the reader stays live between
    components.
    """

    br = NativeBitReader.for_record(packet, rec["payload_offset"], rec["byte_length"])
    rows: list[list[int]] = []
    states: list[dict[str, int]] = []
    selectors = rec["selectors"] or []
    for ci in range(components):
        shift = selectors[ci] if ci < len(selectors) else 0
        for _ in range(row_multiplier):
            if rec["byte_length"] <= 0:
                coeffs = [0] * (groups * 4)
            else:
                coeffs, _widths, _initial = decode_record_component(br, groups, shift)
            rows.append(coeffs)
            states.append({"ptr": br.ptr, "bit": br.bit, "words_left": br.words_left, "status": br.status})
    return rows, states


def infer_packet_width(group: int, packet_type: int, raw_width: int = 7040) -> int:
    """Entropy row width in coefficients for this ARW6 sample's packet groups.

    Type-3 packets carry wavelet detail subbands. Their entropy rows are
    half-width relative to the reconstructed scale reported by the block
    object; the horizontal synthesis stage doubles that later.
    """

    if packet_type == 1:
        return raw_width // 2 if group == 4 else raw_width // 16
    if packet_type == 3:
        if 1 <= group <= 3:
            return raw_width // (1 << (5 - group))
    raise ValueError(f"width inference not yet known for group {group}, type {packet_type}")


def decode_packet_components(
    arw: Path, group: int, index: int, out_prefix: Path | None = None, stream_index: int = 0
) -> dict[str, object]:
    """Decode a type-1 or type-3 packet into one or three int32 component arrays."""

    packet, info = load_packet(arw, group, index, stream_index=stream_index)
    packet_type = info["type2"]
    components = 1 if packet_type == 1 else 3
    width = infer_packet_width(group, packet_type, int(info.get("work_width", info.get("raw_width", 7040))))
    row_multiplier = packet_row_multiplier(group)
    groups_per_row = (width + 3) // 4
    planes: list[list[list[int]]] = [[] for _ in range(components)]
    final_states: list[list[dict[str, int]]] = []
    for rec in info["records"]:
        comp_rows, states = decode_record_components(packet, rec, groups_per_row, components, row_multiplier)
        for ci in range(components):
            start = ci * row_multiplier
            end = start + row_multiplier
            planes[ci].extend(row[:width] for row in comp_rows[start:end])
        final_states.append(states)
    outs: list[str] = []
    if out_prefix:
        import numpy as np

        out_prefix.parent.mkdir(parents=True, exist_ok=True)
        for ci, rows in enumerate(planes):
            arr = np.asarray(rows, dtype=np.int32)
            path = out_prefix.with_name(f"{out_prefix.name}_c{ci}.bin")
            arr.tofile(path)
            outs.append(str(path))
    return {
        "packet": {"group": group, "index": index, "type": packet_type},
        "shape": [len(planes[0]), width],
        "row_multiplier": row_multiplier,
        "components": components,
        "first_row_first16": [plane[0][:16] for plane in planes],
        "last_nonempty_row_first16": [
            next((plane[i][:16] for i in range(len(plane) - 1, -1, -1) if any(plane[i])), []) for plane in planes
        ],
        "final_states_tail": final_states[-5:],
        "outs": outs,
    }


def decode_packet_arrays(arw: Path, group: int, index: int, stream_index: int = 0) -> tuple[list["object"], dict[str, object]]:
    """Decode a packet and return its int32 component arrays in memory."""

    import numpy as np

    packet, info = load_packet(arw, group, index, stream_index=stream_index)
    packet_type = info["type2"]
    components = 1 if packet_type == 1 else 3
    width = infer_packet_width(group, packet_type, int(info.get("work_width", info.get("raw_width", 7040))))
    row_multiplier = packet_row_multiplier(group)
    groups_per_row = (width + 3) // 4
    planes: list[list[list[int]]] = [[] for _ in range(components)]
    final_states: list[list[dict[str, int]]] = []
    for rec in info["records"]:
        comp_rows, states = decode_record_components(packet, rec, groups_per_row, components, row_multiplier)
        for ci in range(components):
            start = ci * row_multiplier
            end = start + row_multiplier
            planes[ci].extend(row[:width] for row in comp_rows[start:end])
        final_states.append(states)

    arrays = [np.asarray(rows, dtype=np.int32) for rows in planes]
    meta = {
        "packet": {"group": group, "index": index, "type": packet_type},
        "stream_index": stream_index,
        "shape": [int(arrays[0].shape[0]), int(arrays[0].shape[1])],
        "row_multiplier": row_multiplier,
        "components": components,
        "final_state": final_states[-1] if final_states else [],
    }
    return arrays, meta


def decode_type1_packet(arw: Path, group: int, index: int, out: Path | None = None, stream_index: int = 0) -> dict[str, object]:
    """Decode an independently parsed type-1 packet into int32 coefficient rows."""

    packet, info = load_packet(arw, group, index, stream_index=stream_index)
    if info["type2"] != 1:
        raise ValueError(f"packet g{group}i{index} is type {info['type2']}, not type 1")
    width = infer_packet_width(group, info["type2"], int(info.get("work_width", info.get("raw_width", 7040))))
    groups_per_row = (width + 3) // 4
    rows: list[list[int]] = []
    final_states: list[dict[str, int]] = []
    for rec in info["records"]:
        if rec["byte_length"] <= 0:
            coeffs = [0] * (groups_per_row * 4)
            br = NativeBitReader.for_record(packet, rec["payload_offset"], rec["byte_length"])
        else:
            coeffs, _widths, br, _initial = decode_type1_row(packet, rec, groups_per_row)
        rows.append(coeffs[:width])
        final_states.append({"ptr": br.ptr, "bit": br.bit, "words_left": br.words_left, "status": br.status})
    if out:
        import numpy as np

        arr = np.asarray(rows, dtype=np.int32)
        out.parent.mkdir(parents=True, exist_ok=True)
        arr.tofile(out)
    return {
        "packet": {"group": group, "index": index, "type": info["type2"]},
        "shape": [len(rows), width],
        "first_row_first16": rows[0][:16],
        "last_nonempty_row_first16": next((rows[i][:16] for i in range(len(rows) - 1, -1, -1) if any(rows[i])), []),
        "final_states_tail": final_states[-5:],
        "out": str(out) if out else "",
    }


def integrate_type1_coefficients(coeffs: "object", dc_offset: int) -> "object":
    """Apply the row postprocess after the entropy call at 0x1a7fb3.

    For each row:

      acc0 = int16(coeff[0]) * 2
      out[0] = acc0 >> 1
      acc_i = acc_{i-1} + int16(coeff[i]) * 2
      out[i] = acc_i >> 1

    Applying this to group0/index0 and adding 2048 reproduces the native v4 c0
    lowpass plane.
    """

    import numpy as np

    c = np.asarray(coeffs, dtype=np.int32)
    out = np.empty_like(c, dtype=np.int32)
    signed = c.astype(np.int16).astype(np.int32)
    acc = signed[:, 0] * 2
    out[:, 0] = acc >> 1
    for x in range(1, c.shape[1]):
        acc = acc + signed[:, x] * 2
        out[:, x] = acc >> 1
    return (out + dc_offset).astype(np.int32)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("arw", nargs="?", default="DSC00089.ARW")
    ap.add_argument("--group", type=int, default=0)
    ap.add_argument("--index", type=int, default=2)
    ap.add_argument("--row", type=int, default=280)
    ap.add_argument("--groups", type=int, default=16)
    ap.add_argument("--packet", action="store_true", help="decode the whole type-1 packet instead of one row")
    ap.add_argument("--components", action="store_true", help="decode a type-1/type-3 packet into component files")
    ap.add_argument("--out", default="")
    ns = ap.parse_args()
    if ns.components:
        result = decode_packet_components(Path(ns.arw), ns.group, ns.index, Path(ns.out) if ns.out else None)
    elif ns.packet:
        result = decode_type1_packet(Path(ns.arw), ns.group, ns.index, Path(ns.out) if ns.out else None)
    else:
        result = replay_row(Path(ns.arw), ns.group, ns.index, ns.row, ns.groups)
    text = json.dumps(result, indent=2)
    if ns.out and not ns.packet and not ns.components:
        Path(ns.out).write_text(text, encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
`,qs=`#!/usr/bin/env python3
"""Integer helpers from the LLVC3 reversing notes.

Closer to a lab notebook than a polished codec module. The small functions
below are named around the traces they came from, then reused by the pure
decoder once a stage lines up.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


INTERNAL_BIAS = 2048


def clamp_sample(x: np.ndarray, bits: int) -> np.ndarray:
    return np.clip(x, 0, (1 << bits) - 1)


def merge_average_detail(lo0: np.ndarray, lo1: np.ndarray, detail: np.ndarray, bits: int = 16) -> np.ndarray:
    """Inverse merge from Edit.exe RVA 0x1a9a40.

    The scalar path is simple enough to leave here as a breadcrumb:

        avg = (lo0 + lo1) >> 1
        sample = avg + 2 * signed_detail
        sample = clamp(sample, 0, (1 << bits) - 1)

    Native output buffers are 16-bit, even when the TIFF tags describe 14-bit
    sensor samples, so the helper returns uint16 too.
    """

    x = ((lo0.astype(np.int32) + lo1.astype(np.int32)) >> 1) + 2 * detail.astype(np.int32)
    return clamp_sample(x, bits).astype(np.uint16)


def add_detail(low: np.ndarray, detail: np.ndarray, bits: int = 16) -> np.ndarray:
    """Lowpass plus signed detail, from the 0x1aafd0 helper family."""

    x = low.astype(np.int32) + detail.astype(np.int32)
    return clamp_sample(x, bits).astype(np.uint16)


def add_double_detail(low: np.ndarray, detail: np.ndarray, bits: int = 16) -> np.ndarray:
    """Lowpass plus 2*detail, seen around the 0x1ab2b0 helpers."""

    x = low.astype(np.int32) + 2 * detail.astype(np.int32)
    return clamp_sample(x, bits).astype(np.uint16)


def sony_inv53_1d(low: np.ndarray, high: np.ndarray, axis: int) -> np.ndarray:
    """One-axis inverse 5/3 lifting in Sony's signed working domain."""

    lo = np.asarray(low, dtype=np.int32)
    hi = np.asarray(high, dtype=np.int32)
    if axis == 0:
        hi_prev = np.vstack([hi[:1], hi[:-1]])
        lo2 = lo - ((hi_prev + hi + 2) >> 2)
        lo_next = np.vstack([lo2[1:], lo2[-1:]])
        hi2 = hi + ((lo2 + lo_next) >> 1)
        out = np.empty((lo.shape[0] * 2, lo.shape[1]), dtype=np.int32)
        out[0::2] = lo2
        out[1::2] = hi2
        return out
    if axis == 1:
        hi_prev = np.concatenate([hi[:, :1], hi[:, :-1]], axis=1)
        lo2 = lo - ((hi_prev + hi + 2) >> 2)
        lo_next = np.concatenate([lo2[:, 1:], lo2[:, -1:]], axis=1)
        hi2 = hi + ((lo2 + lo_next) >> 1)
        out = np.empty((lo.shape[0], lo.shape[1] * 2), dtype=np.int32)
        out[:, 0::2] = lo2
        out[:, 1::2] = hi2
        return out
    raise ValueError("axis must be 0 or 1")


def sony_inv53_1d_high_leading(low: np.ndarray, high: np.ndarray) -> np.ndarray:
    """Vertical inverse 5/3 where Sony's guard line makes the high row lead."""

    lo = np.asarray(low, dtype=np.int32)
    hi = np.asarray(high, dtype=np.int32)
    if hi.shape[1] != lo.shape[1]:
        raise ValueError(f"unexpected high-leading shapes: low={lo.shape}, high={hi.shape}")

    if hi.shape[0] == lo.shape[0] + 1:
        lo2 = lo - ((hi[:-1] + hi[1:] + 2) >> 2)
        hi2 = np.empty_like(hi)
        hi2[0] = hi[0] + lo2[0]
        if lo2.shape[0] > 1:
            hi2[1:-1] = hi[1:-1] + ((lo2[:-1] + lo2[1:]) >> 1)
        hi2[-1] = hi[-1] + lo2[-1]
        out = np.empty((lo.shape[0] * 2 + 1, lo.shape[1]), dtype=np.int32)
    elif hi.shape[0] == lo.shape[0]:
        hi_next = np.vstack([hi[1:], hi[-1:]])
        lo2 = lo - ((hi + hi_next + 2) >> 2)
        hi2 = np.empty_like(hi)
        hi2[0] = hi[0] + lo2[0]
        if lo2.shape[0] > 1:
            hi2[1:] = hi[1:] + ((lo2[:-1] + lo2[1:]) >> 1)
        out = np.empty((lo.shape[0] * 2, lo.shape[1]), dtype=np.int32)
    else:
        raise ValueError(f"unexpected high-leading row counts: low={lo.shape}, high={hi.shape}")

    out[0::2] = hi2
    out[1::2] = lo2
    return out


def llvc3_edge_detail(x: np.ndarray, edge_mode: str = "even") -> np.ndarray:
    """Expand Sony's edge-only HH detail row."""

    xi = x.astype(np.int32)
    signed_half_step = np.where(xi > 0, 1, np.where(xi < 0, -1, 0))
    if edge_mode == "even":
        mask = ((xi & 1) == 0) & (xi != 0)
    elif edge_mode == "odd":
        mask = (xi & 1) != 0
    else:
        raise ValueError(f"unknown LLVC3 edge_mode {edge_mode!r}")
    return (2 * xi + np.where(mask, signed_half_step, 0)).astype(np.int32)


def synthesize_llvc3_level(ll: np.ndarray, sub0: np.ndarray, sub1: np.ndarray, sub2: np.ndarray) -> np.ndarray:
    """Synthesize one LLVC3 scale from LL plus three detail subbands.

    Packet component mapping, as verified against Imaging Edge for this ARW6
    sample:

    * sub0: horizontal detail for the low vertical branch (HL)
    * sub1: vertical detail for the low horizontal branch (LH)
    * sub2: diagonal detail (HH)

    The annoying bit is the extra flush row. At the bottom edge Sony feeds it
    across the branches: LH's final row comes from sub0[-1], HH's final row from
    sub1[-1].
    """

    ll_i = np.asarray(ll, dtype=np.int32)
    h = ll_i.shape[0]
    if ll_i.shape[1] != sub0.shape[1] or sub0.shape != sub1.shape or sub1.shape != sub2.shape:
        raise ValueError(f"unexpected subband shapes: ll={ll_i.shape}, sub0={sub0.shape}, sub1={sub1.shape}, sub2={sub2.shape}")
    if sub0.shape[0] < h + 1:
        raise ValueError(f"subbands need one flush row: ll={ll_i.shape}, sub={sub0.shape}")

    lh = np.empty((h, ll_i.shape[1]), dtype=np.int32)
    lh[:-1] = sub1[1:h]
    lh[-1] = sub0[h]

    hh = np.empty((h, ll_i.shape[1]), dtype=np.int32)
    hh[:-1] = sub2[1:h]
    hh[-1] = sub1[h]

    low_horizontal = sony_inv53_1d(ll_i, lh, axis=0)
    high_horizontal = sony_inv53_1d(sub0[:h], hh, axis=0)
    return sony_inv53_1d(low_horizontal, high_horizontal, axis=1)


def trunc_div2(x: np.ndarray) -> np.ndarray:
    """Integer division by two with C/C++ truncation toward zero."""

    a = np.asarray(x, dtype=np.int32)
    return np.where(a >= 0, a // 2, -((-a) // 2)).astype(np.int32)


def synthesize_llvc3_level_stride(
    ll: np.ndarray,
    sub0: np.ndarray,
    sub1: np.ndarray,
    sub2: np.ndarray,
    edge_rows: int,
    bottom_hh_extra: np.ndarray | None = None,
    edge_mode: str = "even",
) -> np.ndarray:
    """Same synthesis, with the larger line-flush padding used above group 1.

    \`edge_rows\` is the row multiplier minus one for the higher scales:
    group1 -> 0, group2 -> 1, group3 -> 2.  The group1 case falls back to the
    smaller helper above.
    """

    if edge_rows == 0:
        return synthesize_llvc3_level(ll, sub0, sub1, sub2)

    ll_i = np.asarray(ll, dtype=np.int32)
    h, w = ll_i.shape
    if sub0.shape != sub1.shape or sub1.shape != sub2.shape or sub0.shape[1] != w:
        raise ValueError(f"unexpected subband shapes: ll={ll_i.shape}, sub0={sub0.shape}, sub1={sub1.shape}, sub2={sub2.shape}")
    if sub0.shape[0] < h + edge_rows * 2:
        raise ValueError(f"not enough line-flush rows: ll={ll_i.shape}, sub={sub0.shape}, edge_rows={edge_rows}")

    hl = np.empty((h, w), dtype=np.int32)
    hl[:edge_rows] = sub0[:edge_rows]
    hl[edge_rows : h - edge_rows] = sub0[2 * edge_rows : h]
    hl[h - edge_rows :] = sub0[h : h + edge_rows]

    lh = np.empty((h, w), dtype=np.int32)
    lh[:edge_rows] = sub0[edge_rows : 2 * edge_rows]
    lh[edge_rows : h - edge_rows] = sub1[2 * edge_rows : h]
    lh[h - edge_rows :] = sub0[h + edge_rows : h + 2 * edge_rows]

    hh = np.empty((h, w), dtype=np.int32)

    hh[:edge_rows] = llvc3_edge_detail(sub1[:edge_rows], edge_mode)
    hh[edge_rows : h - edge_rows] = sub2[np.arange(edge_rows, h - edge_rows) + edge_rows]
    bottom_hh = llvc3_edge_detail(sub1[h : h + edge_rows], edge_mode)
    if bottom_hh_extra is not None:
        extra = np.asarray(bottom_hh_extra, dtype=np.int32)
        if extra.shape != bottom_hh.shape:
            raise ValueError(f"bottom_hh_extra shape {extra.shape} != bottom edge {bottom_hh.shape}")
        bottom_hh = bottom_hh + extra
    hh[h - edge_rows :] = bottom_hh

    low_horizontal = sony_inv53_1d(ll_i, lh, axis=0)
    high_horizontal = sony_inv53_1d(hl, hh, axis=0)
    return sony_inv53_1d(low_horizontal, high_horizontal, axis=1)


def synthesize_llvc3_guard_group1(ll: np.ndarray, sub0: np.ndarray, sub1: np.ndarray, sub2: np.ndarray) -> np.ndarray:
    """Guard-row group 1 synthesis used by non-16-aligned LLVC heights."""

    ll_i = np.asarray(ll, dtype=np.int32)
    h, w = ll_i.shape
    if sub0.shape != sub1.shape or sub1.shape != sub2.shape or sub0.shape[1] != w:
        raise ValueError(f"unexpected subband shapes: ll={ll_i.shape}, sub0={sub0.shape}, sub1={sub1.shape}, sub2={sub2.shape}")
    if sub0.shape[0] < h + 2:
        raise ValueError(f"not enough guarded group1 rows: ll={ll_i.shape}, sub={sub0.shape}")

    lh = np.empty((h + 1, w), dtype=np.int32)
    lh[:-1] = sub1[1 : 1 + h]
    lh[-1] = sub0[h + 1]

    hh = np.empty((h + 1, w), dtype=np.int32)
    hh[:-1] = sub2[1 : 1 + h]
    hh[-1] = sub1[h + 1]

    low_horizontal = sony_inv53_1d_high_leading(ll_i, lh)
    high_horizontal = sony_inv53_1d_high_leading(sub0[1 : 1 + h], hh)
    return sony_inv53_1d(low_horizontal, high_horizontal, axis=1)


def synthesize_llvc3_guard_group2(
    ll: np.ndarray, sub0: np.ndarray, sub1: np.ndarray, sub2: np.ndarray, edge_mode: str = "even"
) -> np.ndarray:
    """Guard-row group 2 synthesis for cropped-height ARW6 tiles."""

    ll_i = np.asarray(ll, dtype=np.int32)
    h, w = ll_i.shape
    if sub0.shape != sub1.shape or sub1.shape != sub2.shape or sub0.shape[1] != w:
        raise ValueError(f"unexpected subband shapes: ll={ll_i.shape}, sub0={sub0.shape}, sub1={sub1.shape}, sub2={sub2.shape}")
    if sub0.shape[0] < h + 3:
        raise ValueError(f"not enough guarded group2 rows: ll={ll_i.shape}, sub={sub0.shape}")

    hl = sub0[2 : 2 + h]
    lh = np.empty((h, w), dtype=np.int32)
    lh[0] = sub0[0]
    lh[1:] = sub1[2 : 1 + h]

    hh = np.empty((h, w), dtype=np.int32)
    hh[0] = llvc3_edge_detail(sub0[1:2], edge_mode)[0]
    hh[1:] = sub2[2 : 1 + h]

    low_horizontal = sony_inv53_1d_high_leading(ll_i, lh)
    high_horizontal = sony_inv53_1d_high_leading(hl, hh)
    return sony_inv53_1d(low_horizontal, high_horizontal, axis=1)


def synthesize_llvc3_guard_group3(
    ll: np.ndarray, sub0: np.ndarray, sub1: np.ndarray, sub2: np.ndarray, edge_mode: str = "even"
) -> np.ndarray:
    """Guard-row group 3 synthesis for cropped-height ARW6 tiles."""

    ll_i = np.asarray(ll, dtype=np.int32)
    h, w = ll_i.shape
    if sub0.shape != sub1.shape or sub1.shape != sub2.shape or sub0.shape[1] != w:
        raise ValueError(f"unexpected subband shapes: ll={ll_i.shape}, sub0={sub0.shape}, sub1={sub1.shape}, sub2={sub2.shape}")
    if sub0.shape[0] < h + 5:
        raise ValueError(f"not enough guarded group3 rows: ll={ll_i.shape}, sub={sub0.shape}")

    hl = np.empty((h, w), dtype=np.int32)
    hl[0] = sub0[0]
    hl[1:] = sub0[4 : 3 + h]

    lh = np.empty((h, w), dtype=np.int32)
    lh[0] = sub0[1]
    lh[1:-1] = sub1[4 : 2 + h]
    lh[-1] = sub0[h + 3]

    hh = np.empty((h, w), dtype=np.int32)
    hh[0] = llvc3_edge_detail(sub0[2:3], edge_mode)[0]
    hh[1:-1] = sub2[4 : 2 + h]
    hh[-1] = llvc3_edge_detail(sub0[h + 4 : h + 5], edge_mode)[0]

    low_horizontal = sony_inv53_1d(ll_i, lh, axis=0)
    high_horizontal = sony_inv53_1d(hl, hh, axis=0)
    return sony_inv53_1d(low_horizontal, high_horizontal, axis=1)


def synthesize_llvc3_final_green(ll: np.ndarray, detail: np.ndarray, top_rows: int = 4) -> np.ndarray:
    """Final CFA-green reconstruction, from the 0x1ab570 path.

    Not the same 2-D 5/3 inverse used by groups 1..3. It expands the half-width
    green lowpass into both RGGB green sites. The row offsets look odd because
    the native line buffer keeps four guard rows at the top and a few latency
    rows at the bottom.
    """

    ll_i = np.asarray(ll, dtype=np.int32)
    det = np.asarray(detail, dtype=np.int32)
    h, w = ll_i.shape
    if not 0 <= top_rows <= 8:
        raise ValueError(f"unexpected final green top row count {top_rows}")
    if det.shape[1] != w or det.shape[0] < 8 + max(0, h - top_rows):
        raise ValueError(f"unexpected final green shapes: ll={ll_i.shape}, detail={det.shape}")

    selected = np.empty((h, w), dtype=np.int32)
    top = min(top_rows, h)
    selected[:top] = det[:top]
    if h > top:
        selected[top:] = det[8 : 8 + (h - top)]

    odd_green = np.empty((h, w), dtype=np.int32)
    for y in range(h):
        cur = selected[y]
        prev = selected[y - 1] if y > 0 else cur
        pred = np.empty(w, dtype=np.int32)
        pred[:-1] = (cur[1:] + prev[:-1] + cur[:-1] + prev[1:]) >> 2
        pred[-1] = ((prev[-1] + cur[-1]) * 2) >> 2
        odd_green[y] = ((2 * ll_i[y] - pred) >> 1).astype(np.int32)

    even_green = np.empty((h, w), dtype=np.int32)
    for y in range(h):
        cur = odd_green[y]
        nxt = odd_green[y + 1] if y + 1 < h else cur
        even_green[y, 0] = selected[y, 0] + (((cur[0] + nxt[0]) * 2) >> 2)
        if w > 1:
            even_green[y, 1:] = selected[y, 1:] + ((cur[:-1] + nxt[:-1] + nxt[1:] + cur[1:]) >> 2)

    out = np.empty((h, w * 2), dtype=np.int32)
    out[:, 0::2] = even_green
    out[:, 1::2] = odd_green
    return out


def finalize_llvc3_color_planes(
    v1_green: np.ndarray, v1_red: np.ndarray, v1_blue: np.ndarray, full_green: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Apply the final color-plane relation for the decoded output.

    Red and blue keep the group-3 residuals. Only the green predictor is swapped
    for the final CFA-green average:

        residual = (v1_color - v1_green) / 2
        v0_color = average(clamp12(final_green_pair)) + 2 * residual

    Sony clamps the two final green sites to the 12-bit code domain before using
    them as the red/blue predictor. Without that, highlight green overshoot leaks
    into R/B as one to three LUT code steps.
    """

    g = np.asarray(v1_green, dtype=np.int32)
    r = np.asarray(v1_red, dtype=np.int32)
    b = np.asarray(v1_blue, dtype=np.int32)
    fg = np.asarray(full_green, dtype=np.int32)
    if fg.shape != (g.shape[0], g.shape[1] * 2):
        raise ValueError(f"unexpected final green shape: v1={g.shape}, full={fg.shape}")
    if r.shape != g.shape or b.shape != g.shape:
        raise ValueError(f"unexpected v1 color shapes: green={g.shape}, red={r.shape}, blue={b.shape}")

    fg_pred = np.clip(fg + INTERNAL_BIAS, 0, 4095).astype(np.int32) - INTERNAL_BIAS
    gavg = (fg_pred[:, 0::2] + fg_pred[:, 1::2]) >> 1
    red_residual = (r - g) // 2
    blue_residual = (b - g) // 2
    return fg, gavg + 2 * red_residual, gavg + 2 * blue_residual


def signed_to_sample(x: np.ndarray, bits: int = 16, bias: int = INTERNAL_BIAS) -> np.ndarray:
    """Convert LLVC3 signed internal rows to Sony's unsigned output samples."""

    return clamp_sample(np.asarray(x, dtype=np.int32) + bias, bits).astype(np.uint16)


def apply_sample_lut(code_samples: np.ndarray, lut: np.ndarray) -> np.ndarray:
    """Map unsigned LLVC3 code-domain samples through a Sony sample LUT."""

    table = np.asarray(lut, dtype=np.uint16).reshape(-1)
    if table.size == 0:
        raise ValueError("sample LUT is empty")
    code = np.clip(np.asarray(code_samples, dtype=np.int32), 0, table.size - 1)
    return table[code].astype(np.uint16)


def clamp_signed_to_code_range(x: np.ndarray, max_code: int = 4095, bias: int = INTERNAL_BIAS) -> np.ndarray:
    """Clamp signed LLVC3 rows to Sony's 12-bit code range, then return signed rows."""

    return np.clip(np.asarray(x, dtype=np.int32) + bias, 0, max_code).astype(np.int32) - bias


def lifting_predict_detail(detail: np.ndarray, a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray) -> np.ndarray:
    """Prediction/update kernel from the 0x1ab570 scalar path.

    The four neighbor names are still placeholder-ish, but the integer
    operation itself is clear in the trace:

        pred = (a + b + c + d) >> 2
        out = (2 * detail - pred) >> 1
    """

    pred = (a.astype(np.int32) + b.astype(np.int32) + c.astype(np.int32) + d.astype(np.int32)) >> 2
    return ((2 * detail.astype(np.int32) - pred) >> 1).astype(np.int32)


def recombine_rggb(c0: np.ndarray, c1: np.ndarray, c2: np.ndarray) -> np.ndarray:
    """Recombine decoded LLVC planes into the TIFF-declared RGGB Bayer mosaic."""

    half_h, width = c0.shape
    if c1.shape != (half_h, width // 2) or c2.shape != (half_h, width // 2):
        raise ValueError(f"unexpected plane shapes: {c0.shape}, {c1.shape}, {c2.shape}")
    out = np.empty((half_h * 2, width), dtype=np.uint16)
    out[0::2, 0::2] = c1
    out[0::2, 1::2] = c0[:, 1::2]
    out[1::2, 0::2] = c0[:, 0::2]
    out[1::2, 1::2] = c2
    return out


def selftest() -> dict[str, object]:
    lo0 = np.array([1000, 1002, 10, 65530], dtype=np.uint16)
    lo1 = np.array([1002, 1004, 10, 65530], dtype=np.uint16)
    detail = np.array([0, 1, -20, 20], dtype=np.int32)
    merged = merge_average_detail(lo0, lo1, detail, bits=16)
    added = add_detail(lo0, detail, bits=16)
    doubled = add_double_detail(lo0, detail, bits=16)
    pred = lifting_predict_detail(
        detail,
        np.array([4, 8, 12, 16]),
        np.array([4, 8, 12, 16]),
        np.array([4, 8, 12, 16]),
        np.array([4, 8, 12, 16]),
    )
    return {
        "merge_average_detail": merged.tolist(),
        "add_detail": added.tolist(),
        "add_double_detail": doubled.tolist(),
        "lifting_predict_detail": pred.tolist(),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--out", default="")
    ns = ap.parse_args()
    result = selftest()
    text = json.dumps(result, indent=2)
    if ns.out:
        Path(ns.out).write_text(text, encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
`,Ks="/assets/joraw-BAR7UfD7.js",$s="/assets/joraw-CZg3hVRg.wasm";const yr=512,Js=Ws["Sony ILCE-7RM5"].split(/\s+/).map(Number).filter(Number.isFinite);let ki=null,Fi=null;function Zs(n){if(n.byteLength<8)return null;const t=n.getUint16(0,!1);return t===18761?!0:t===19789?!1:null}function Kr(n,t,e){const i=Math.min(n.length,t+e);let r="";for(let s=t;s<i;s++){const o=n[s];if(o===0)break;r+=String.fromCharCode(o)}return r.trim()}function jt(n,t,e,i,r){const s=e===1||e===2||e===7?1:e===3||e===8?2:e===4||e===9?4:0;if(!s)return[];const o=s*i,a=o<=4?r:n.getUint32(r,t);if(a<0||a+o>n.byteLength)return[];const l=[];for(let u=0;u<i;u++){const h=a+u*s;e===1||e===2||e===7?l.push(n.getUint8(h)):e===3?l.push(n.getUint16(h,t)):e===8?l.push(n.getInt16(h,t)):e===4?l.push(n.getUint32(h,t)):e===9&&l.push(n.getInt32(h,t))}return l}function xr(n,t,e,i,r,s){if(i!==2||r<=0)return"";const o=r<=4?s:t.getUint32(s,e);return o<0||o>=n.length?"":Kr(n,o,r)}function $r(n){const t=new Uint8Array(n),e=new DataView(n),i=Zs(e);if(i===null||e.getUint16(2,i)!==42)return null;const s=c=>e.getUint16(c,i),o=c=>e.getUint32(c,i),a=[o(4)],l=new Set;let u="",h="";for(;a.length;){const c=a.pop();if(l.has(c)||c<=0||c+2>e.byteLength)continue;l.add(c);const m=s(c);if(c+2+m*12+4>e.byteLength)continue;const d=new Map;for(let _=0;_<m;_++){const P=c+2+_*12,S=s(P),w=s(P+2),C=o(P+4),A=P+8;d.set(S,{type:w,count:C,valueOffset:A})}const f=d.get(271),p=d.get(272);f&&!u&&(u=xr(t,e,i,f.type,f.count,f.valueOffset)),p&&!h&&(h=xr(t,e,i,p.type,p.count,p.valueOffset));const g=d.get(330);if(g){const _=jt(e,i,g.type,g.count,g.valueOffset);for(const P of _)a.push(P)}const y=d.get(259),x=d.get(262);if(y&&x){const _=jt(e,i,y.type,y.count,y.valueOffset)[0],P=jt(e,i,x.type,x.count,x.valueOffset)[0];if(_===32766&&P===32803){const S=jt(e,i,d.get(256).type,d.get(256).count,d.get(256).valueOffset)[0],w=jt(e,i,d.get(257).type,d.get(257).count,d.get(257).valueOffset)[0],C=jt(e,i,d.get(258).type,d.get(258).count,d.get(258).valueOffset)[0],A=jt(e,i,d.get(273).type,d.get(273).count,d.get(273).valueOffset)[0],k=jt(e,i,d.get(279).type,d.get(279).count,d.get(279).valueOffset)[0],F=d.get(33422)?jt(e,i,d.get(33422).type,d.get(33422).count,d.get(33422).valueOffset):[0,1,1,2];d.get(29456)&&jt(e,i,d.get(29456).type,d.get(29456).count,d.get(29456).valueOffset);const M=d.get(50717)?jt(e,i,d.get(50717).type,d.get(50717).count,d.get(50717).valueOffset)[0]:16383,v=d.get(50719)?jt(e,i,d.get(50719).type,d.get(50719).count,d.get(50719).valueOffset):[],I=d.get(50720)?jt(e,i,d.get(50720).type,d.get(50720).count,d.get(50720).valueOffset):[];if(A+yr+16>t.length||A+k>t.length)return null;const R=A+yr,U=Kr(t,R,4),Y=t[R+8]<<8|t[R+9],B=t[R+10]<<8|t[R+11],L=t[R+12]<<8|t[R+13],O=t[R+14]<<8|t[R+15],D=L>>4&63,X=O>>13,et=O>>10&3,V=B*2,$=k>=4?(t[A]|t[A+1]<<8|t[A+2]<<16|t[A+3]<<24)>>>0:0,H=Y===S&&V===w;let it=!1;if($>=1&&$<=16&&k>=8+$*24){const Q=new Map,T=new Map;let j=!0;for(let G=0;G<$;G++){const z=A+8+G*24,Z=e.getUint32(z+8,!0),tt=e.getUint32(z+12,!0),nt=e.getUint32(z+16,!0),ct=e.getUint32(z+20,!0);if(!nt||!ct||Z+nt>S||tt+ct>w){j=!1;break}const dt=Q.get(tt);if(dt!==void 0&&dt!==ct){j=!1;break}Q.set(tt,ct),T.set(tt,(T.get(tt)||0)+nt)}if(j){const G=Array.from(Q.keys()).sort((Z,tt)=>Z-tt);let z=0;for(const Z of G){if(Z!==z||T.get(Z)!==S){j=!1;break}z+=Q.get(Z)}it=j&&z===w}}const at=$>=1&&$<=16&&Y>0&&V>0&&S%Y===0&&w%V===0&&$===S/Y*(w/V);if(U!=="A000"&&U!=="0000"||!H&&!at&&!it||D!==16||X!==3||et!==3)return null;const E=[1024,1024,1024,1024],J=F.slice(0,4).map(Q=>Q===0?"R":Q===2?"B":"G").join("")||"RGGB";return{width:S,height:w,bitsPerSample:C,compression:_,photometric:P,blackLevel:E,whiteLevel:Number(M||16383),cfaPattern:J,defaultCropOrigin:v.length>=2?[Number(v[0]),Number(v[1])]:void 0,defaultCropSize:I.length>=2?[Number(I[0]),Number(I[1])]:void 0,make:u||"SONY",model:h||"ILCE-7M5"}}}const b=o(c+2+m*12);b&&a.push(b)}return null}async function ta(n){return ki||(ki=(async()=>{const t=await Hs();return t.__jtrSonyCrawHqDecoderReady||(await t.FS.mkdirTree("/sony_craw_hq"),await t.FS.writeFile("/sony_craw_hq/llvc3_bitstream_probe.py",js),await t.FS.writeFile("/sony_craw_hq/llvc3_entropy.py",Qs),await t.FS.writeFile("/sony_craw_hq/llvc3_math.py",qs),await t.runPythonAsync(`
import sys
if "/sony_craw_hq" not in sys.path:
    sys.path.insert(0, "/sony_craw_hq")
from pathlib import Path
import numpy as np
from llvc3_bitstream_probe import find_llvc_streams, find_raw_subifd
from llvc3_entropy import decode_packet_arrays, integrate_type1_coefficients
from llvc3_math import apply_sample_lut, clamp_signed_to_code_range, recombine_rggb, signed_to_sample
from llvc3_math import finalize_llvc3_color_planes, synthesize_llvc3_final_green, synthesize_llvc3_level_stride
from llvc3_math import synthesize_llvc3_guard_group1, synthesize_llvc3_guard_group2, synthesize_llvc3_guard_group3

def jtr_align_up(value, multiple):
    return ((value + multiple - 1) // multiple) * multiple

def jtr_validate_stream_layout(raw_info, streams):
    if not streams:
        raise ValueError("no LLVC3 streams to validate")
    rows = {}
    for stream in streams:
        if int(stream.tile_x) < 0 or int(stream.tile_y) < 0:
            raise ValueError(f"negative LLVC3 tile position: {stream}")
        if int(stream.tile_width) != int(stream.header.coded_width):
            raise ValueError(f"LLVC3 tile width/header mismatch: {stream}")
        if int(stream.tile_height) != int(stream.header.logical_height):
            raise ValueError(f"LLVC3 tile height/header mismatch: {stream}")
        if int(stream.tile_x) + int(stream.tile_width) > int(raw_info.width):
            raise ValueError(f"LLVC3 tile exceeds raw width: {stream}")
        if int(stream.tile_y) + int(stream.tile_height) > int(raw_info.height):
            raise ValueError(f"LLVC3 tile exceeds raw height: {stream}")
        rows.setdefault(int(stream.tile_y), []).append(stream)

    expected_y = 0
    for tile_y in sorted(rows):
        row = sorted(rows[tile_y], key=lambda s: int(s.tile_x))
        row_height = int(row[0].tile_height)
        if tile_y != expected_y:
            raise ValueError(f"LLVC3 tile rows have a gap/overlap at y={tile_y}, expected {expected_y}")
        expected_x = 0
        for stream in row:
            if int(stream.tile_height) != row_height:
                raise ValueError(f"LLVC3 row has mixed tile heights at y={tile_y}")
            if int(stream.tile_x) != expected_x:
                raise ValueError(f"LLVC3 tile columns have a gap/overlap at x={stream.tile_x}, expected {expected_x}")
            expected_x += int(stream.tile_width)
        if expected_x != int(raw_info.width):
            raise ValueError(f"LLVC3 tile row width {expected_x} does not cover raw width {raw_info.width}")
        expected_y += row_height
    if expected_y != int(raw_info.height):
        raise ValueError(f"LLVC3 tile rows height {expected_y} does not cover raw height {raw_info.height}")

def jtr_combine_tiled_arrays(tiles, streams, x_divisor=1, fill=0):
    if len(tiles) != len(streams):
        raise ValueError(f"{len(tiles)} decoded tiles do not match {len(streams)} LLVC3 streams")
    ys = sorted({int(stream.tile_y) for stream in streams})
    y_rank = {y: i for i, y in enumerate(ys)}
    placements = []
    max_x = 0
    max_y = 0
    for tile, stream in zip(tiles, streams):
        x = int(stream.tile_x) // x_divisor
        y = y_rank[int(stream.tile_y)] * tile.shape[0]
        placements.append((x, y, tile))
        max_x = max(max_x, x + tile.shape[1])
        max_y = max(max_y, y + tile.shape[0])
    out = np.full((max_y, max_x), fill, dtype=tiles[0].dtype)
    for x, y, tile in placements:
        out[y : y + tile.shape[0], x : x + tile.shape[1]] = tile
    return out

def jtr_decode_signed_planes(arw, stream_index=0, stream_header=None):
    if stream_header is None:
        raw_info, strip = find_raw_subifd(arw)
        streams = find_llvc_streams(strip)
        if not streams:
            raise ValueError("no LLVC3 stream found in ARW6 raw strip")
        stream_header = streams[stream_index].header
    coded_height = stream_header.logical_height
    padded_height = jtr_align_up(coded_height, 16)
    guarded_height = coded_height != padded_height
    low_rows = padded_height // 16
    low_start = 1 if guarded_height else 0
    low_count = low_rows - low_start

    g0, _meta = decode_packet_arrays(arw, 0, 0, stream_index=stream_index)
    green = integrate_type1_coefficients(g0[0][low_start : low_start + low_count], 2048) - 2048

    r0, _meta = decode_packet_arrays(arw, 0, 1, stream_index=stream_index)
    red_residual = integrate_type1_coefficients(r0[0][low_start : low_start + low_count], 0)

    b0, _meta = decode_packet_arrays(arw, 0, 2, stream_index=stream_index)
    blue_residual = integrate_type1_coefficients(b0[0][low_start : low_start + low_count], 0)

    for group, edge_rows in ((1, 0), (2, 1), (3, 2)):
        old_green = green
        old_red_residual = red_residual
        old_blue_residual = blue_residual

        planes, _meta = decode_packet_arrays(arw, group, 0, stream_index=stream_index)
        if guarded_height:
            if group == 1:
                green = synthesize_llvc3_guard_group1(old_green, planes[0], planes[1], planes[2])
            elif group == 2:
                green = synthesize_llvc3_guard_group2(old_green, planes[0], planes[1], planes[2])
            else:
                green = synthesize_llvc3_guard_group3(old_green, planes[0], planes[1], planes[2])
        else:
            green = synthesize_llvc3_level_stride(old_green, planes[0], planes[1], planes[2], edge_rows)

        planes, _meta = decode_packet_arrays(arw, group, 1, stream_index=stream_index)
        edge_mode = "odd" if group == 3 else "even"
        if guarded_height:
            if group == 1:
                red_residual = synthesize_llvc3_guard_group1(old_red_residual, planes[0], planes[1], planes[2])
            elif group == 2:
                red_residual = synthesize_llvc3_guard_group2(old_red_residual, planes[0], planes[1], planes[2])
            else:
                red_residual = synthesize_llvc3_guard_group3(
                    old_red_residual, planes[0], planes[1], planes[2], edge_mode=edge_mode
                )
        else:
            red_residual = synthesize_llvc3_level_stride(
                old_red_residual, planes[0], planes[1], planes[2], edge_rows, edge_mode=edge_mode
            )

        planes, _meta = decode_packet_arrays(arw, group, 2, stream_index=stream_index)
        if guarded_height:
            if group == 1:
                blue_residual = synthesize_llvc3_guard_group1(old_blue_residual, planes[0], planes[1], planes[2])
            elif group == 2:
                blue_residual = synthesize_llvc3_guard_group2(old_blue_residual, planes[0], planes[1], planes[2])
            else:
                blue_residual = synthesize_llvc3_guard_group3(
                    old_blue_residual, planes[0], planes[1], planes[2], edge_mode=edge_mode
                )
        else:
            blue_residual = synthesize_llvc3_level_stride(
                old_blue_residual, planes[0], planes[1], planes[2], edge_rows, edge_mode=edge_mode
            )

    g4, _meta = decode_packet_arrays(arw, 4, 0, stream_index=stream_index)
    full_green = synthesize_llvc3_final_green(green, g4[0], top_rows=2 if guarded_height else 4)
    v1_red = green + 2 * red_residual
    v1_blue = green + 2 * blue_residual
    c0, c1, c2 = finalize_llvc3_color_planes(green, v1_red, v1_blue, full_green)
    if c0.shape[0] != stream_header.coded_half_height:
        extra_rows = c0.shape[0] - stream_header.coded_half_height
        if extra_rows < 0:
            raise ValueError(
                f"stream {stream_index} decoded only {c0.shape[0]} half-height rows, "
                f"expected {stream_header.coded_half_height}"
            )
        crop_top = 0 if guarded_height else extra_rows // 2
        bottom = crop_top + stream_header.coded_half_height
        c0 = c0[crop_top:bottom]
        c1 = c1[crop_top:bottom]
        c2 = c2[crop_top:bottom]
    return c0, c1, c2

def jtr_decode_sony_craw_hq(arw_bytes, lut_bytes=None):
    path = Path("/tmp/jtr_sony_craw_hq_input.arw")
    path.write_bytes(bytes(arw_bytes))
    raw_info, strip = find_raw_subifd(path)
    streams = find_llvc_streams(strip)
    if not streams:
        raise ValueError("no LLVC3 stream found in ARW6 raw strip")
    if any(s.header.component_count != 3 for s in streams):
        raise ValueError(f"unexpected ARW6/LLVC3 component count in streams: {streams}")
    jtr_validate_stream_layout(raw_info, streams)
    if raw_info.width % 16 or raw_info.height % 16:
        raise ValueError(f"decoder expects dimensions divisible by 16, got {raw_info.width}x{raw_info.height}")

    signed_tiles = []
    for stream_index, stream in enumerate(streams):
        signed_c0, signed_c1, signed_c2 = jtr_decode_signed_planes(path, stream_index, stream.header)
        if signed_c0.shape != (stream.header.coded_half_height, stream.header.coded_width):
            raise ValueError(
                f"stream {stream_index} c0 decoded to {signed_c0.shape}, expected "
                f"{stream.header.coded_half_height}x{stream.header.coded_width}"
            )
        expected_chroma_shape = (stream.header.coded_half_height, stream.header.coded_width // 2)
        if signed_c1.shape != expected_chroma_shape or signed_c2.shape != expected_chroma_shape:
            raise ValueError(
                f"stream {stream_index} chroma decoded to {signed_c1.shape}/{signed_c2.shape}, "
                f"expected {expected_chroma_shape}"
            )
        signed_tiles.append((signed_c0, signed_c1, signed_c2))

    if lut_bytes is not None:
        lut = np.frombuffer(bytes(lut_bytes), dtype="<u2")
        if lut.size:
            if lut.size < 65536:
                lut = np.pad(lut, (0, 65536 - lut.size), constant_values=int(lut[-1]))
            lut = lut[:65536].astype(np.uint16)
        else:
            lut = None
    else:
        lut = None

    tile_raws = []
    for signed_c0, signed_c1, signed_c2 in signed_tiles:
        if lut is not None:
            sample_c0 = apply_sample_lut(signed_to_sample(clamp_signed_to_code_range(signed_c0)), lut)
            sample_c1 = apply_sample_lut(signed_to_sample(clamp_signed_to_code_range(signed_c1)), lut)
            sample_c2 = apply_sample_lut(signed_to_sample(clamp_signed_to_code_range(signed_c2)), lut)
        else:
            sample_c0 = signed_to_sample(signed_c0)
            sample_c1 = signed_to_sample(signed_c1)
            sample_c2 = signed_to_sample(signed_c2)
        tile_raws.append(recombine_rggb(sample_c0, sample_c1, sample_c2))
    raw = jtr_combine_tiled_arrays(tile_raws, streams, x_divisor=1, fill=1024)
    if raw.shape != (raw_info.height, raw_info.width):
        raise ValueError(f"decoded raw shape {raw.shape} does not match TIFF raw {raw_info.height}x{raw_info.width}")
    return raw.astype("<u2", copy=False).tobytes()
`),t.__jtrSonyCrawHqDecoderReady=!0),t})()),ki}function ea(n){return $r(n)}async function na(){return Fi||(Fi=import(Ks).then(n=>{const t=n.default;if(typeof t!="function")throw new Error("Sony cRAW HQ LibRaw WASM module is missing its initializer");return t({locateFile:(e,i)=>e.endsWith("joraw.wasm")?$s:i+e})})),Fi}function ia(n){return n==="ILCE-7M5"?Js:null}function Jr(n,t,e,i){var l;if(n.length!==t.width*t.height)throw new Error(`Sony cRAW HQ decoded size mismatch: got ${n.length}, expected ${t.width*t.height}`);const r=t.model||"ILCE-7M5",s=r.startsWith("Sony ")?r:`Sony ${r}`,o=ia(r),a={...i||{},make:t.make||(i==null?void 0:i.camera_make)||"SONY",model:r,camera_make:t.make||(i==null?void 0:i.camera_make)||"SONY",camera_model:r,UniqueCameraModel:s,sourceFormat:e==="libraw-wasm"?"Sony cRAW HQ / LLVC3 (LibRaw WASM)":"Sony cRAW HQ / LLVC3",sonyCrawHq:{...t,decodeBackend:e},color_desc:t.cfaPattern,black_level_per_channel:t.blackLevel,white_level:t.whiteLevel,color_matrix:o&&o.length===9?o:void 0,idata:{filters:2492765332,colors:3},color_data:{...(i==null?void 0:i.color_data)||{},black:1024,cblack_rawpy_style:t.blackLevel,dng_levels:{...((l=i==null?void 0:i.color_data)==null?void 0:l.dng_levels)||{},dng_cblack:t.blackLevel,dng_whitelevel:t.whiteLevel}}};return{data:n,width:t.width,height:t.height,bayerPattern:t.cfaPattern,blackLevels:t.blackLevel,whiteLevel:t.whiteLevel,metadata:a,isThreePlane:!1,isXTrans:!1}}async function ra(n,t,e){const i=typeof performance<"u"?performance.now():Date.now(),r=await na(),s=typeof performance<"u"?performance.now():Date.now(),o=r.LibRaw||r.JoRaw;if(!o)throw new Error("Sony cRAW HQ LibRaw WASM class not found");const a=new o;try{const l=new Uint8Array(n);a.open(l,{});const u=typeof performance<"u"?performance.now():Date.now();let h=null;try{h=a.metadata(!0)}catch(y){console.warn("[Sony cRAW HQ] fast WASM metadata read failed",y)}const c=typeof performance<"u"?performance.now():Date.now(),m=a.getRawImage(),d=typeof performance<"u"?performance.now():Date.now();if(!m||!m.data)throw new Error("Sony cRAW HQ LibRaw WASM returned no raw image");const f=m.data instanceof Uint16Array?m.data:new Uint16Array(m.data.buffer,m.data.byteOffset||0,m.data.byteLength/2),p=typeof performance<"u"?performance.now():Date.now();if(m.width!==t.width||m.height!==t.height)throw new Error(`Sony cRAW HQ LibRaw WASM dimensions mismatch: got ${m.width}x${m.height}, expected ${t.width}x${t.height}`);const g=typeof performance<"u"?performance.now():Date.now();return console.info("[Sony cRAW HQ] fast decode timings",{width:t.width,height:t.height,backend:"libraw-wasm",wasmReadyMs:Math.round(s-i),openMs:Math.round(u-s),metadataMs:Math.round(c-u),unpackMs:Math.round(d-c),copyMs:Math.round(p-d),totalMs:Math.round(g-i)}),{rawImageData:Jr(f,t,"libraw-wasm",h),info:t}}finally{typeof a.delete=="function"?a.delete():typeof a.close=="function"&&a.close()}}async function sa(n,t,e){const i=typeof performance<"u"?performance.now():Date.now(),r=await ta(),s=typeof performance<"u"?performance.now():Date.now(),o=new Uint8Array(n),a=await fetch(new URL("/assets/sony_llvc3_static_lut4096_padded_u16-FsVBk-IV.bin",import.meta.url));if(!a.ok)throw new Error(`Failed to load Sony LLVC3 sample LUT: HTTP ${a.status}`);const l=new Uint8Array(await a.arrayBuffer()),u=typeof performance<"u"?performance.now():Date.now();r.globals.set("jtr_sony_arw_bytes",o),r.globals.set("jtr_sony_lut_bytes",l);const h=await r.runPythonAsync("jtr_decode_sony_craw_hq(jtr_sony_arw_bytes.to_py(), jtr_sony_lut_bytes.to_py())"),c=typeof performance<"u"?performance.now():Date.now(),m=h.toJs();typeof h.destroy=="function"&&h.destroy(),r.globals.delete("jtr_sony_arw_bytes"),r.globals.delete("jtr_sony_lut_bytes");const d=new Uint8Array(m.byteLength);d.set(m);const f=new Uint16Array(d.buffer),p=typeof performance<"u"?performance.now():Date.now(),g=typeof performance<"u"?performance.now():Date.now();return console.info("[Sony cRAW HQ] decode timings",{width:t.width,height:t.height,backend:"pyodide",pyodideReadyMs:Math.round(s-i),lutLoadMs:Math.round(u-s),llvc3DecodeMs:Math.round(c-u),copyMs:Math.round(p-c),totalMs:Math.round(g-i)}),{rawImageData:Jr(f,t,"pyodide"),info:t}}async function aa(n,t){const e=$r(n);if(!e)return null;try{return await ra(n,e,t)}catch(i){return console.warn("[Sony cRAW HQ] fast WASM decode failed; falling back to Pyodide",i),sa(n,e)}}async function oa(n,t){return aa(n,t)}const la={1:1,2:1,3:2,4:4,5:8,7:1,9:4,10:8,11:4,12:8},ca=34713,Oi=155;let Ai=null;function br(n,t,e){return un(n,t,e,Math.min(e.count,256)).map(r=>r&255).filter(r=>r!==0).map(r=>String.fromCharCode(r)).join("").trim()}function un(n,t,e,i=Number.MAX_SAFE_INTEGER){const r=la[e.type]||0;if(!r)return[];const o=e.count*r<=4?e.valueOffset:n.getUint32(e.valueOffset,t),a=Math.min(e.count,i),l=[];for(let u=0;u<a;u++){const h=o+u*r;if(h<0||h+r>n.byteLength)break;switch(e.type){case 1:case 2:case 7:l.push(n.getUint8(h));break;case 3:l.push(n.getUint16(h,t));break;case 4:l.push(n.getUint32(h,t));break;case 9:l.push(n.getInt32(h,t));break;case 5:{const c=n.getUint32(h,t),m=n.getUint32(h+4,t);l.push(m?c/m:0);break}case 10:{const c=n.getInt32(h,t),m=n.getInt32(h+4,t);l.push(m?c/m:0);break}}}return l}function _r(n,t,e){if(e<=0||e+2>n.byteLength)return null;const i=n.getUint16(e,t),r=2+i*12+4;if(i>4096||e+r>n.byteLength)return null;const s=new Map;for(let o=0;o<i;o++){const a=e+2+o*12,l=n.getUint16(a,t);s.set(l,{tag:l,type:n.getUint16(a+2,t),count:n.getUint32(a+4,t),valueOffset:a+8})}return s}function ua(n,t,e){if(e<=0||e+2>n.byteLength)return 0;const i=n.getUint16(e,t),r=e+2+i*12;return r+4>n.byteLength?0:n.getUint32(r,t)}function ha(n){const t=["R","G","B","C","M","Y","W"];return n.length<4?"RGGB":n.slice(0,4).map(e=>t[e]||"G").join("")}function Qe(n,t,e,i){const r=n.get(i);if(!r)return null;const s=un(t,e,r,1);return s.length?s[0]:null}function Zr(n){if(n.byteLength<16)return null;const t=new DataView(n),e=t.getUint16(0,!1),i=e===18761;if(!i&&e!==19789||t.getUint16(2,i)!==42)return null;const r=t.getUint32(4,i),s=_r(t,i,r);if(!s)return null;const o=s.get(271),a=s.get(272),l=o?br(t,i,o):"",u=a?br(t,i,a):"";if(!/NIKON/i.test(l+" "+u))return null;const h=Qe(s,t,i,274)||1,c=s.get(330);if(!c)return null;const m=un(t,i,c,16);for(const f of m){const p=_r(t,i,f);if(!p)continue;const g=Qe(p,t,i,259)||0,y=Qe(p,t,i,273)||0,x=Qe(p,t,i,279)||0,b=Qe(p,t,i,256)||0,_=Qe(p,t,i,257)||0;if(g!==ca||!y||!x||b<=0||_<=0||y+x>n.byteLength)continue;const P=Qe(p,t,i,258)||14,S=p.get(33422),w=S?ha(un(t,i,S,4)):"RGGB",C=p.get(50714),A=p.get(50717),k=C?un(t,i,C,4):[],F=A?un(t,i,A,1):[],M=k.length&&Number(k[0])||0,v=[Number(k[0])||M,Number(k[1])||M,Number(k[2])||M,Number(k[3])||M],I=y+Oi,R=I+3<n.byteLength?t.getUint8(I+3):0;return{make:l,model:u,width:b,height:_,stripOffset:y,stripByteCount:x,bayerPattern:w,blackLevels:v,whiteLevel:Math.max(1,Number(F[0])||(1<<Math.min(P,15))-1||16383),orientation:h,compression:g,firstBp:R,isHeStar:R>=1&&R<=2}}return ua(t,i,r)>0,null}function da(){return{wasi_snapshot_preview1:{fd_close:()=>0,fd_seek:(n,t,e,i,r)=>0,fd_write:()=>0,proc_exit:n=>{throw new Error(`Nikon HE WASM exited with code ${n}`)}}}}async function fa(){return Ai||(Ai=(async()=>{const n=new URL("/assets/nikon-he-decoder-DUHNoGN2.wasm",import.meta.url),t=await fetch(n);if(!t.ok)throw new Error(`Failed to load Nikon HE WASM: ${t.status}`);const e=await t.arrayBuffer(),r=(await WebAssembly.instantiate(e,da())).instance.exports;if(!r.memory||!r.nikon_he_decode||!r.nikon_he_malloc||!r.nikon_he_free)throw new Error("Nikon HE WASM exports are incomplete.");return r})()),Ai}function pa(n){return Zr(n)}async function ma(n){const t=Zr(n);if(!t)return null;const e=t.stripOffset+Oi,i=t.stripByteCount-Oi;if(i<=12||e+i>n.byteLength)throw new Error("Invalid Nikon HE precinct bounds.");const r=await fa(),s=t.width*t.height*2,o=r.nikon_he_malloc(i),a=r.nikon_he_malloc(s);try{new Uint8Array(r.memory.buffer).set(new Uint8Array(n,e,i),o);const u=r.nikon_he_decode(o,i,t.width,t.height,a);if(u<0)throw new Error(`Nikon HE decoder failed with code ${u}.`);const h=new Uint16Array(r.memory.buffer,a,t.width*t.height),c=new Uint16Array(h);let m=65535,d=0;const f=Math.max(1,Math.floor(c.length/2e5));for(let p=0;p<c.length;p+=f){const g=c[p];g<m&&(m=g),g>d&&(d=g)}if(m===d)throw new Error("Nikon HE decoder returned a constant image; this HE variant is not supported by the current WASM core.");return{data:c,width:t.width,height:t.height,bayerPattern:t.bayerPattern,blackLevels:t.blackLevels,whiteLevel:t.whiteLevel,metadata:{description:`Nikon ${t.isHeStar?"HE*":"HE"} RAW decoded by standalone WASM`,camera_make:t.make,camera_model:t.model,make:t.make,model:t.model,orientation:t.orientation,compression:t.compression,raw_width:t.width,raw_height:t.height,width:t.width,height:t.height,color_desc:t.bayerPattern,white_level:t.whiteLevel,black_level_per_channel:t.blackLevels,nikon_he:{isHeStar:t.isHeStar,firstBp:t.firstBp,stripOffset:t.stripOffset,stripByteCount:t.stripByteCount,precinctOffset:e,precinctSize:i}},isThreePlane:!1}}finally{r.nikon_he_free(a),r.nikon_he_free(o)}}var ri=typeof self<"u"?self:global;const kn=typeof navigator<"u",ga=kn&&typeof HTMLImageElement>"u",Kn=!(typeof global>"u"||typeof process>"u"||!process.versions||!process.versions.node),si=ri.Buffer,Vn=ri.BigInt,ai=!!si,ya=n=>n;function $n(n,t=ya){if(Kn)try{return typeof require=="function"?Promise.resolve(t(require(n))):import(n).then(t)}catch{console.warn(`Couldn't load ${n}`)}}let Qi=ri.fetch;const xa=n=>Qi=n;if(!ri.fetch){const n=$n("http",i=>i),t=$n("https",i=>i),e=(i,{headers:r}={})=>new Promise(async(s,o)=>{let{port:a,hostname:l,pathname:u,protocol:h,search:c}=new URL(i);const m={method:"GET",hostname:l,path:encodeURI(u)+c,headers:r};a!==""&&(m.port=Number(a));const d=(h==="https:"?await t:await n).request(m,f=>{if(f.statusCode===301||f.statusCode===302){let p=new URL(f.headers.location,i).toString();return e(p,{headers:r}).then(s).catch(o)}s({status:f.statusCode,arrayBuffer:()=>new Promise(p=>{let g=[];f.on("data",y=>g.push(y)),f.on("end",()=>p(Buffer.concat(g)))})})});d.on("error",o),d.end()});xa(e)}function st(n,t,e){return t in n?Object.defineProperty(n,t,{value:e,enumerable:!0,configurable:!0,writable:!0}):n[t]=e,n}const Jn=n=>ts(n)?void 0:n,ba=n=>n!==void 0;function ts(n){return n===void 0||(n instanceof Map?n.size===0:Object.values(n).filter(ba).length===0)}function kt(n){let t=new Error(n);throw delete t.stack,t}function Ke(n){return(n=function(t){for(;t.endsWith("\0");)t=t.slice(0,-1);return t}(n).trim())===""?void 0:n}function Gi(n){let t=function(e){let i=0;return e.ifd0.enabled&&(i+=1024),e.exif.enabled&&(i+=2048),e.makerNote&&(i+=2048),e.userComment&&(i+=1024),e.gps.enabled&&(i+=512),e.interop.enabled&&(i+=100),e.ifd1.enabled&&(i+=1024),i+2048}(n);return n.jfif.enabled&&(t+=50),n.xmp.enabled&&(t+=2e4),n.iptc.enabled&&(t+=14e3),n.icc.enabled&&(t+=6e3),t}const Vi=n=>String.fromCharCode.apply(null,n),Mr=typeof TextDecoder<"u"?new TextDecoder("utf-8"):void 0;function es(n){return Mr?Mr.decode(n):ai?Buffer.from(n).toString("utf8"):decodeURIComponent(escape(Vi(n)))}class Wt{static from(t,e){return t instanceof this&&t.le===e?t:new Wt(t,void 0,void 0,e)}constructor(t,e=0,i,r){if(typeof r=="boolean"&&(this.le=r),Array.isArray(t)&&(t=new Uint8Array(t)),t===0)this.byteOffset=0,this.byteLength=0;else if(t instanceof ArrayBuffer){i===void 0&&(i=t.byteLength-e);let s=new DataView(t,e,i);this._swapDataView(s)}else if(t instanceof Uint8Array||t instanceof DataView||t instanceof Wt){i===void 0&&(i=t.byteLength-e),(e+=t.byteOffset)+i>t.byteOffset+t.byteLength&&kt("Creating view outside of available memory in ArrayBuffer");let s=new DataView(t.buffer,e,i);this._swapDataView(s)}else if(typeof t=="number"){let s=new DataView(new ArrayBuffer(t));this._swapDataView(s)}else kt("Invalid input argument for BufferView: "+t)}_swapArrayBuffer(t){this._swapDataView(new DataView(t))}_swapBuffer(t){this._swapDataView(new DataView(t.buffer,t.byteOffset,t.byteLength))}_swapDataView(t){this.dataView=t,this.buffer=t.buffer,this.byteOffset=t.byteOffset,this.byteLength=t.byteLength}_lengthToEnd(t){return this.byteLength-t}set(t,e,i=Wt){return t instanceof DataView||t instanceof Wt?t=new Uint8Array(t.buffer,t.byteOffset,t.byteLength):t instanceof ArrayBuffer&&(t=new Uint8Array(t)),t instanceof Uint8Array||kt("BufferView.set(): Invalid data argument."),this.toUint8().set(t,e),new i(this,e,t.byteLength)}subarray(t,e){return e=e||this._lengthToEnd(t),new Wt(this,t,e)}toUint8(){return new Uint8Array(this.buffer,this.byteOffset,this.byteLength)}getUint8Array(t,e){return new Uint8Array(this.buffer,this.byteOffset+t,e)}getString(t=0,e=this.byteLength){return es(this.getUint8Array(t,e))}getLatin1String(t=0,e=this.byteLength){let i=this.getUint8Array(t,e);return Vi(i)}getUnicodeString(t=0,e=this.byteLength){const i=[];for(let r=0;r<e&&t+r<this.byteLength;r+=2)i.push(this.getUint16(t+r));return Vi(i)}getInt8(t){return this.dataView.getInt8(t)}getUint8(t){return this.dataView.getUint8(t)}getInt16(t,e=this.le){return this.dataView.getInt16(t,e)}getInt32(t,e=this.le){return this.dataView.getInt32(t,e)}getUint16(t,e=this.le){return this.dataView.getUint16(t,e)}getUint32(t,e=this.le){return this.dataView.getUint32(t,e)}getFloat32(t,e=this.le){return this.dataView.getFloat32(t,e)}getFloat64(t,e=this.le){return this.dataView.getFloat64(t,e)}getFloat(t,e=this.le){return this.dataView.getFloat32(t,e)}getDouble(t,e=this.le){return this.dataView.getFloat64(t,e)}getUintBytes(t,e,i){switch(e){case 1:return this.getUint8(t,i);case 2:return this.getUint16(t,i);case 4:return this.getUint32(t,i);case 8:return this.getUint64&&this.getUint64(t,i)}}getUint(t,e,i){switch(e){case 8:return this.getUint8(t,i);case 16:return this.getUint16(t,i);case 32:return this.getUint32(t,i);case 64:return this.getUint64&&this.getUint64(t,i)}}toString(t){return this.dataView.toString(t,this.constructor.name)}ensureChunk(){}}function Xi(n,t){kt(`${n} '${t}' was not loaded, try using full build of exifr.`)}class qi extends Map{constructor(t){super(),this.kind=t}get(t,e){return this.has(t)||Xi(this.kind,t),e&&(t in e||function(i,r){kt(`Unknown ${i} '${r}'.`)}(this.kind,t),e[t].enabled||Xi(this.kind,t)),super.get(t)}keyList(){return Array.from(this.keys())}}var ye=new qi("file parser"),Pt=new qi("segment parser"),be=new qi("file reader");function _a(n,t){return typeof n=="string"?wr(n,t):kn&&!ga&&n instanceof HTMLImageElement?wr(n.src,t):n instanceof Uint8Array||n instanceof ArrayBuffer||n instanceof DataView?new Wt(n):kn&&n instanceof Blob?zi(n,t,"blob",sn):void kt("Invalid input argument")}function wr(n,t){return(e=n).startsWith("data:")||e.length>1e4?Yi(n,t,"base64"):Kn&&n.includes("://")?zi(n,t,"url",rn):Kn?Yi(n,t,"fs"):kn?zi(n,t,"url",rn):void kt("Invalid input argument");var e}async function zi(n,t,e,i){return be.has(e)?Yi(n,t,e):i?async function(r,s){let o=await s(r);return new Wt(o)}(n,i):void kt(`Parser ${e} is not loaded`)}async function Yi(n,t,e){let i=new(be.get(e))(n,t);return await i.read(),i}const rn=n=>Qi(n).then(t=>t.arrayBuffer()),sn=n=>new Promise((t,e)=>{let i=new FileReader;i.onloadend=()=>t(i.result||new ArrayBuffer),i.onerror=e,i.readAsArrayBuffer(n)});class Ma extends Map{get tagKeys(){return this.allKeys||(this.allKeys=Array.from(this.keys())),this.allKeys}get tagValues(){return this.allValues||(this.allValues=Array.from(this.values())),this.allValues}}function Ct(n,t,e){let i=new Ma;for(let[r,s]of e)i.set(r,s);if(Array.isArray(t))for(let r of t)n.set(r,i);else n.set(t,i);return i}function an(n,t,e){let i,r=n.get(t);for(i of e)r.set(i[0],i[1])}const Nt=new Map,oe=new Map,Ee=new Map,Te=["chunked","firstChunkSize","firstChunkSizeNode","firstChunkSizeBrowser","chunkSize","chunkLimit"],gn=["jfif","xmp","icc","iptc","ihdr"],on=["tiff",...gn],St=["ifd0","ifd1","exif","gps","interop"],Ie=[...on,...St],Ne=["makerNote","userComment"],yn=["translateKeys","translateValues","reviveValues","multiSegment"],Re=[...yn,"sanitize","mergeOutput","silentErrors"];class ns{get translate(){return this.translateKeys||this.translateValues||this.reviveValues}}class Sn extends ns{get needed(){return this.enabled||this.deps.size>0}constructor(t,e,i,r){if(super(),st(this,"enabled",!1),st(this,"skip",new Set),st(this,"pick",new Set),st(this,"deps",new Set),st(this,"translateKeys",!1),st(this,"translateValues",!1),st(this,"reviveValues",!1),this.key=t,this.enabled=e,this.parse=this.enabled,this.applyInheritables(r),this.canBeFiltered=St.includes(t),this.canBeFiltered&&(this.dict=Nt.get(t)),i!==void 0)if(Array.isArray(i))this.parse=this.enabled=!0,this.canBeFiltered&&i.length>0&&this.translateTagSet(i,this.pick);else if(typeof i=="object"){if(this.enabled=!0,this.parse=i.parse!==!1,this.canBeFiltered){let{pick:s,skip:o}=i;s&&s.length>0&&this.translateTagSet(s,this.pick),o&&o.length>0&&this.translateTagSet(o,this.skip)}this.applyInheritables(i)}else i===!0||i===!1?this.parse=this.enabled=i:kt(`Invalid options argument: ${i}`)}applyInheritables(t){let e,i;for(e of yn)i=t[e],i!==void 0&&(this[e]=i)}translateTagSet(t,e){if(this.dict){let i,r,{tagKeys:s,tagValues:o}=this.dict;for(i of t)typeof i=="string"?(r=o.indexOf(i),r===-1&&(r=s.indexOf(Number(i))),r!==-1&&e.add(Number(s[r]))):e.add(i)}else for(let i of t)e.add(i)}finalizeFilters(){!this.enabled&&this.deps.size>0?(this.enabled=!0,Zn(this.pick,this.deps)):this.enabled&&this.pick.size>0&&Zn(this.pick,this.deps)}}var Ot={jfif:!1,tiff:!0,xmp:!1,icc:!1,iptc:!1,ifd0:!0,ifd1:!1,exif:!0,gps:!0,interop:!1,ihdr:void 0,makerNote:!1,userComment:!1,multiSegment:!1,skip:[],pick:[],translateKeys:!0,translateValues:!0,reviveValues:!0,sanitize:!0,mergeOutput:!0,silentErrors:!0,chunked:!0,firstChunkSize:void 0,firstChunkSizeNode:512,firstChunkSizeBrowser:65536,chunkSize:65536,chunkLimit:5},Sr=new Map;class ln extends ns{static useCached(t){let e=Sr.get(t);return e!==void 0||(e=new this(t),Sr.set(t,e)),e}constructor(t){super(),t===!0?this.setupFromTrue():t===void 0?this.setupFromUndefined():Array.isArray(t)?this.setupFromArray(t):typeof t=="object"?this.setupFromObject(t):kt(`Invalid options argument ${t}`),this.firstChunkSize===void 0&&(this.firstChunkSize=kn?this.firstChunkSizeBrowser:this.firstChunkSizeNode),this.mergeOutput&&(this.ifd1.enabled=!1),this.filterNestedSegmentTags(),this.traverseTiffDependencyTree(),this.checkLoadedPlugins()}setupFromUndefined(){let t;for(t of Te)this[t]=Ot[t];for(t of Re)this[t]=Ot[t];for(t of Ne)this[t]=Ot[t];for(t of Ie)this[t]=new Sn(t,Ot[t],void 0,this)}setupFromTrue(){let t;for(t of Te)this[t]=Ot[t];for(t of Re)this[t]=Ot[t];for(t of Ne)this[t]=!0;for(t of Ie)this[t]=new Sn(t,!0,void 0,this)}setupFromArray(t){let e;for(e of Te)this[e]=Ot[e];for(e of Re)this[e]=Ot[e];for(e of Ne)this[e]=Ot[e];for(e of Ie)this[e]=new Sn(e,!1,void 0,this);this.setupGlobalFilters(t,void 0,St)}setupFromObject(t){let e;for(e of(St.ifd0=St.ifd0||St.image,St.ifd1=St.ifd1||St.thumbnail,Object.assign(this,t),Te))this[e]=Ti(t[e],Ot[e]);for(e of Re)this[e]=Ti(t[e],Ot[e]);for(e of Ne)this[e]=Ti(t[e],Ot[e]);for(e of on)this[e]=new Sn(e,Ot[e],t[e],this);for(e of St)this[e]=new Sn(e,Ot[e],t[e],this.tiff);this.setupGlobalFilters(t.pick,t.skip,St,Ie),t.tiff===!0?this.batchEnableWithBool(St,!0):t.tiff===!1?this.batchEnableWithUserValue(St,t):Array.isArray(t.tiff)?this.setupGlobalFilters(t.tiff,void 0,St):typeof t.tiff=="object"&&this.setupGlobalFilters(t.tiff.pick,t.tiff.skip,St)}batchEnableWithBool(t,e){for(let i of t)this[i].enabled=e}batchEnableWithUserValue(t,e){for(let i of t){let r=e[i];this[i].enabled=r!==!1&&r!==void 0}}setupGlobalFilters(t,e,i,r=i){if(t&&t.length){for(let o of r)this[o].enabled=!1;let s=vr(t,i);for(let[o,a]of s)Zn(this[o].pick,a),this[o].enabled=!0}else if(e&&e.length){let s=vr(e,i);for(let[o,a]of s)Zn(this[o].skip,a)}}filterNestedSegmentTags(){let{ifd0:t,exif:e,xmp:i,iptc:r,icc:s}=this;this.makerNote?e.deps.add(37500):e.skip.add(37500),this.userComment?e.deps.add(37510):e.skip.add(37510),i.enabled||t.skip.add(700),r.enabled||t.skip.add(33723),s.enabled||t.skip.add(34675)}traverseTiffDependencyTree(){let{ifd0:t,exif:e,gps:i,interop:r}=this;r.needed&&(e.deps.add(40965),t.deps.add(40965)),e.needed&&t.deps.add(34665),i.needed&&t.deps.add(34853),this.tiff.enabled=St.some(s=>this[s].enabled===!0)||this.makerNote||this.userComment;for(let s of St)this[s].finalizeFilters()}get onlyTiff(){return!gn.map(t=>this[t].enabled).some(t=>t===!0)&&this.tiff.enabled}checkLoadedPlugins(){for(let t of on)this[t].enabled&&!Pt.has(t)&&Xi("segment parser",t)}}function vr(n,t){let e,i,r,s,o=[];for(r of t){for(s of(e=Nt.get(r),i=[],e))(n.includes(s[0])||n.includes(s[1]))&&i.push(s[0]);i.length&&o.push([r,i])}return o}function Ti(n,t){return n!==void 0?n:t!==void 0?t:void 0}function Zn(n,t){for(let e of t)n.add(e)}st(ln,"default",Ot);class Ue{constructor(t){st(this,"parsers",{}),st(this,"output",{}),st(this,"errors",[]),st(this,"pushToErrors",e=>this.errors.push(e)),this.options=ln.useCached(t)}async read(t){this.file=await _a(t,this.options)}setup(){if(this.fileParser)return;let{file:t}=this,e=t.getUint16(0);for(let[i,r]of ye)if(r.canHandle(t,e))return this.fileParser=new r(this.options,this.file,this.parsers),t[i]=!0;this.file.close&&this.file.close(),kt("Unknown file format")}async parse(){let{output:t,errors:e}=this;return this.setup(),this.options.silentErrors?(await this.executeParsers().catch(this.pushToErrors),e.push(...this.fileParser.errors)):await this.executeParsers(),this.file.close&&this.file.close(),this.options.silentErrors&&e.length>0&&(t.errors=e),Jn(t)}async executeParsers(){let{output:t}=this;await this.fileParser.parse();let e=Object.values(this.parsers).map(async i=>{let r=await i.parse();i.assignToOutput(t,r)});this.options.silentErrors&&(e=e.map(i=>i.catch(this.pushToErrors))),await Promise.all(e)}async extractThumbnail(){this.setup();let{options:t,file:e}=this,i=Pt.get("tiff",t);var r;if(e.tiff?r={start:0,type:"tiff"}:e.jpeg&&(r=await this.fileParser.getOrFindSegment("tiff")),r===void 0)return;let s=await this.fileParser.ensureSegmentChunk(r),o=this.parsers.tiff=new i(s,t,e),a=await o.extractThumbnail();return e.close&&e.close(),a}}async function oi(n,t){let e=new Ue(t);return await e.read(n),e.parse()}var wa=Object.freeze({__proto__:null,parse:oi,Exifr:Ue,fileParsers:ye,segmentParsers:Pt,fileReaders:be,tagKeys:Nt,tagValues:oe,tagRevivers:Ee,createDictionary:Ct,extendDictionary:an,fetchUrlAsArrayBuffer:rn,readBlobAsArrayBuffer:sn,chunkedProps:Te,otherSegments:gn,segments:on,tiffBlocks:St,segmentsAndBlocks:Ie,tiffExtractables:Ne,inheritables:yn,allFormatters:Re,Options:ln});class li{constructor(t,e,i){st(this,"errors",[]),st(this,"ensureSegmentChunk",async r=>{let s=r.start,o=r.size||65536;if(this.file.chunked)if(this.file.available(s,o))r.chunk=this.file.subarray(s,o);else try{r.chunk=await this.file.readChunk(s,o)}catch(a){kt(`Couldn't read segment: ${JSON.stringify(r)}. ${a.message}`)}else this.file.byteLength>s+o?r.chunk=this.file.subarray(s,o):r.size===void 0?r.chunk=this.file.subarray(s):kt("Segment unreachable: "+JSON.stringify(r));return r.chunk}),this.extendOptions&&this.extendOptions(t),this.options=t,this.file=e,this.parsers=i}injectSegment(t,e){this.options[t].enabled&&this.createParser(t,e)}createParser(t,e){let i=new(Pt.get(t))(e,this.options,this.file);return this.parsers[t]=i}createParsers(t){for(let e of t){let{type:i,chunk:r}=e,s=this.options[i];if(s&&s.enabled){let o=this.parsers[i];o&&o.append||o||this.createParser(i,r)}}}async readSegments(t){let e=t.map(this.ensureSegmentChunk);await Promise.all(e)}}class ae{static findPosition(t,e){let i=t.getUint16(e+2)+2,r=typeof this.headerLength=="function"?this.headerLength(t,e,i):this.headerLength,s=e+r,o=i-r;return{offset:e,length:i,headerLength:r,start:s,size:o,end:s+o}}static parse(t,e={}){return new this(t,new ln({[this.type]:e}),t).parse()}normalizeInput(t){return t instanceof Wt?t:new Wt(t)}constructor(t,e={},i){st(this,"errors",[]),st(this,"raw",new Map),st(this,"handleError",r=>{if(!this.options.silentErrors)throw r;this.errors.push(r.message)}),this.chunk=this.normalizeInput(t),this.file=i,this.type=this.constructor.type,this.globalOptions=this.options=e,this.localOptions=e[this.type],this.canTranslate=this.localOptions&&this.localOptions.translate}translate(){this.canTranslate&&(this.translated=this.translateBlock(this.raw,this.type))}get output(){return this.translated?this.translated:this.raw?Object.fromEntries(this.raw):void 0}translateBlock(t,e){let i=Ee.get(e),r=oe.get(e),s=Nt.get(e),o=this.options[e],a=o.reviveValues&&!!i,l=o.translateValues&&!!r,u=o.translateKeys&&!!s,h={};for(let[c,m]of t)a&&i.has(c)?m=i.get(c)(m):l&&r.has(c)&&(m=this.translateValue(m,r.get(c))),u&&s.has(c)&&(c=s.get(c)||c),h[c]=m;return h}translateValue(t,e){return e[t]||e.DEFAULT||t}assignToOutput(t,e){this.assignObjectToOutput(t,this.constructor.type,e)}assignObjectToOutput(t,e,i){if(this.globalOptions.mergeOutput)return Object.assign(t,i);t[e]?Object.assign(t[e],i):t[e]=i}}st(ae,"headerLength",4),st(ae,"type",void 0),st(ae,"multiSegment",!1),st(ae,"canHandle",()=>!1);function Sa(n){return n===192||n===194||n===196||n===219||n===221||n===218||n===254}function va(n){return n>=224&&n<=239}function Ca(n,t,e){for(let[i,r]of Pt)if(r.canHandle(n,t,e))return i}class Cr extends li{constructor(...t){super(...t),st(this,"appSegments",[]),st(this,"jpegSegments",[]),st(this,"unknownSegments",[])}static canHandle(t,e){return e===65496}async parse(){await this.findAppSegments(),await this.readSegments(this.appSegments),this.mergeMultiSegments(),this.createParsers(this.mergedAppSegments||this.appSegments)}setupSegmentFinderArgs(t){t===!0?(this.findAll=!0,this.wanted=new Set(Pt.keyList())):(t=t===void 0?Pt.keyList().filter(e=>this.options[e].enabled):t.filter(e=>this.options[e].enabled&&Pt.has(e)),this.findAll=!1,this.remaining=new Set(t),this.wanted=new Set(t)),this.unfinishedMultiSegment=!1}async findAppSegments(t=0,e){this.setupSegmentFinderArgs(e);let{file:i,findAll:r,wanted:s,remaining:o}=this;if(!r&&this.file.chunked&&(r=Array.from(s).some(a=>{let l=Pt.get(a),u=this.options[a];return l.multiSegment&&u.multiSegment}),r&&await this.file.readWhole()),t=this.findAppSegmentsInRange(t,i.byteLength),!this.options.onlyTiff&&i.chunked){let a=!1;for(;o.size>0&&!a&&(i.canReadNextChunk||this.unfinishedMultiSegment);){let{nextChunkOffset:l}=i,u=this.appSegments.some(h=>!this.file.available(h.offset||h.start,h.length||h.size));if(a=t>l&&!u?!await i.readNextChunk(t):!await i.readNextChunk(l),(t=this.findAppSegmentsInRange(t,i.byteLength))===void 0)return}}}findAppSegmentsInRange(t,e){e-=2;let i,r,s,o,a,l,{file:u,findAll:h,wanted:c,remaining:m,options:d}=this;for(;t<e;t++)if(u.getUint8(t)===255){if(i=u.getUint8(t+1),va(i)){if(r=u.getUint16(t+2),s=Ca(u,t,r),s&&c.has(s)&&(o=Pt.get(s),a=o.findPosition(u,t),l=d[s],a.type=s,this.appSegments.push(a),!h&&(o.multiSegment&&l.multiSegment?(this.unfinishedMultiSegment=a.chunkNumber<a.chunkCount,this.unfinishedMultiSegment||m.delete(s)):m.delete(s),m.size===0)))break;d.recordUnknownSegments&&(a=ae.findPosition(u,t),a.marker=i,this.unknownSegments.push(a)),t+=r+1}else if(Sa(i)){if(r=u.getUint16(t+2),i===218&&d.stopAfterSos!==!1)return;d.recordJpegSegments&&this.jpegSegments.push({offset:t,length:r,marker:i}),t+=r+1}}return t}mergeMultiSegments(){if(!this.appSegments.some(e=>e.multiSegment))return;let t=function(e,i){let r,s,o,a=new Map;for(let l=0;l<e.length;l++)r=e[l],s=r[i],a.has(s)?o=a.get(s):a.set(s,o=[]),o.push(r);return Array.from(a)}(this.appSegments,"type");this.mergedAppSegments=t.map(([e,i])=>{let r=Pt.get(e,this.options);return r.handleMultiSegments?{type:e,chunk:r.handleMultiSegments(i)}:i[0]})}getSegment(t){return this.appSegments.find(e=>e.type===t)}async getOrFindSegment(t){let e=this.getSegment(t);return e===void 0&&(await this.findAppSegments(0,[t]),e=this.getSegment(t)),e}}st(Cr,"type","jpeg"),ye.set("jpeg",Cr);const Pa=[void 0,1,1,2,4,8,1,1,2,4,8,4,8,4];class ka extends ae{parseHeader(){var t=this.chunk.getUint16();t===18761?this.le=!0:t===19789&&(this.le=!1),this.chunk.le=this.le,this.headerParsed=!0}parseTags(t,e,i=new Map){let{pick:r,skip:s}=this.options[e];r=new Set(r);let o=r.size>0,a=s.size===0,l=this.chunk.getUint16(t);t+=2;for(let u=0;u<l;u++){let h=this.chunk.getUint16(t);if(o){if(r.has(h)&&(i.set(h,this.parseTag(t,h,e)),r.delete(h),r.size===0))break}else!a&&s.has(h)||i.set(h,this.parseTag(t,h,e));t+=12}return i}parseTag(t,e,i){let{chunk:r}=this,s=r.getUint16(t+2),o=r.getUint32(t+4),a=Pa[s];if(a*o<=4?t+=8:t=r.getUint32(t+8),(s<1||s>13)&&kt(`Invalid TIFF value type. block: ${i.toUpperCase()}, tag: ${e.toString(16)}, type: ${s}, offset ${t}`),t>r.byteLength&&kt(`Invalid TIFF value offset. block: ${i.toUpperCase()}, tag: ${e.toString(16)}, type: ${s}, offset ${t} is outside of chunk size ${r.byteLength}`),s===1)return r.getUint8Array(t,o);if(s===2)return Ke(r.getString(t,o));if(s===7)return r.getUint8Array(t,o);if(o===1)return this.parseTagValue(s,t);{let l=new(function(h){switch(h){case 1:return Uint8Array;case 3:return Uint16Array;case 4:return Uint32Array;case 5:return Array;case 6:return Int8Array;case 8:return Int16Array;case 9:return Int32Array;case 10:return Array;case 11:return Float32Array;case 12:return Float64Array;default:return Array}}(s))(o),u=a;for(let h=0;h<o;h++)l[h]=this.parseTagValue(s,t),t+=u;return l}}parseTagValue(t,e){let{chunk:i}=this;switch(t){case 1:return i.getUint8(e);case 3:return i.getUint16(e);case 4:return i.getUint32(e);case 5:return i.getUint32(e)/i.getUint32(e+4);case 6:return i.getInt8(e);case 8:return i.getInt16(e);case 9:return i.getInt32(e);case 10:return i.getInt32(e)/i.getInt32(e+4);case 11:return i.getFloat(e);case 12:return i.getDouble(e);case 13:return i.getUint32(e);default:kt(`Invalid tiff type ${t}`)}}}class Ii extends ka{static canHandle(t,e){return t.getUint8(e+1)===225&&t.getUint32(e+4)===1165519206&&t.getUint16(e+8)===0}async parse(){this.parseHeader();let{options:t}=this;return t.ifd0.enabled&&await this.parseIfd0Block(),t.exif.enabled&&await this.safeParse("parseExifBlock"),t.gps.enabled&&await this.safeParse("parseGpsBlock"),t.interop.enabled&&await this.safeParse("parseInteropBlock"),t.ifd1.enabled&&await this.safeParse("parseThumbnailBlock"),this.createOutput()}safeParse(t){let e=this[t]();return e.catch!==void 0&&(e=e.catch(this.handleError)),e}findIfd0Offset(){this.ifd0Offset===void 0&&(this.ifd0Offset=this.chunk.getUint32(4))}findIfd1Offset(){if(this.ifd1Offset===void 0){this.findIfd0Offset();let t=this.chunk.getUint16(this.ifd0Offset),e=this.ifd0Offset+2+12*t;this.ifd1Offset=this.chunk.getUint32(e)}}parseBlock(t,e){let i=new Map;return this[e]=i,this.parseTags(t,e,i),i}async parseIfd0Block(){if(this.ifd0)return;let{file:t}=this;this.findIfd0Offset(),this.ifd0Offset<8&&kt("Malformed EXIF data"),!t.chunked&&this.ifd0Offset>t.byteLength&&kt(`IFD0 offset points to outside of file.
this.ifd0Offset: ${this.ifd0Offset}, file.byteLength: ${t.byteLength}`),t.tiff&&await t.ensureChunk(this.ifd0Offset,Gi(this.options));let e=this.parseBlock(this.ifd0Offset,"ifd0");return e.size!==0?(this.exifOffset=e.get(34665),this.interopOffset=e.get(40965),this.gpsOffset=e.get(34853),this.xmp=e.get(700),this.iptc=e.get(33723),this.icc=e.get(34675),this.options.sanitize&&(e.delete(34665),e.delete(40965),e.delete(34853),e.delete(700),e.delete(33723),e.delete(34675)),e):void 0}async parseExifBlock(){if(this.exif||(this.ifd0||await this.parseIfd0Block(),this.exifOffset===void 0))return;this.file.tiff&&await this.file.ensureChunk(this.exifOffset,Gi(this.options));let t=this.parseBlock(this.exifOffset,"exif");return this.interopOffset||(this.interopOffset=t.get(40965)),this.makerNote=t.get(37500),this.userComment=t.get(37510),this.options.sanitize&&(t.delete(40965),t.delete(37500),t.delete(37510)),this.unpack(t,41728),this.unpack(t,41729),t}unpack(t,e){let i=t.get(e);i&&i.length===1&&t.set(e,i[0])}async parseGpsBlock(){if(this.gps||(this.ifd0||await this.parseIfd0Block(),this.gpsOffset===void 0))return;let t=this.parseBlock(this.gpsOffset,"gps");return t&&t.has(2)&&t.has(4)&&(t.set("latitude",Pr(...t.get(2),t.get(1))),t.set("longitude",Pr(...t.get(4),t.get(3)))),t}async parseInteropBlock(){if(!this.interop&&(this.ifd0||await this.parseIfd0Block(),this.interopOffset!==void 0||this.exif||await this.parseExifBlock(),this.interopOffset!==void 0))return this.parseBlock(this.interopOffset,"interop")}async parseThumbnailBlock(t=!1){if(!this.ifd1&&!this.ifd1Parsed&&(!this.options.mergeOutput||t))return this.findIfd1Offset(),this.ifd1Offset>0&&(this.parseBlock(this.ifd1Offset,"ifd1"),this.ifd1Parsed=!0),this.ifd1}async extractThumbnail(){if(this.headerParsed||this.parseHeader(),this.ifd1Parsed||await this.parseThumbnailBlock(!0),this.ifd1===void 0)return;let t=this.ifd1.get(513),e=this.ifd1.get(514);return this.chunk.getUint8Array(t,e)}get image(){return this.ifd0}get thumbnail(){return this.ifd1}createOutput(){let t,e,i,r={};for(e of St)if(t=this[e],!ts(t))if(i=this.canTranslate?this.translateBlock(t,e):Object.fromEntries(t),this.options.mergeOutput){if(e==="ifd1")continue;Object.assign(r,i)}else r[e]=i;return this.makerNote&&(r.makerNote=this.makerNote),this.userComment&&(r.userComment=this.userComment),r}assignToOutput(t,e){if(this.globalOptions.mergeOutput)Object.assign(t,e);else for(let[i,r]of Object.entries(e))this.assignObjectToOutput(t,i,r)}}function Pr(n,t,e,i){var r=n+t/60+e/3600;return i!=="S"&&i!=="W"||(r*=-1),r}st(Ii,"type","tiff"),st(Ii,"headerLength",10),Pt.set("tiff",Ii);var Fa=Object.freeze({__proto__:null,default:wa,Exifr:Ue,fileParsers:ye,segmentParsers:Pt,fileReaders:be,tagKeys:Nt,tagValues:oe,tagRevivers:Ee,createDictionary:Ct,extendDictionary:an,fetchUrlAsArrayBuffer:rn,readBlobAsArrayBuffer:sn,chunkedProps:Te,otherSegments:gn,segments:on,tiffBlocks:St,segmentsAndBlocks:Ie,tiffExtractables:Ne,inheritables:yn,allFormatters:Re,Options:ln,parse:oi});const Ki={ifd0:!1,ifd1:!1,exif:!1,gps:!1,interop:!1,sanitize:!1,reviveValues:!0,translateKeys:!1,translateValues:!1,mergeOutput:!1},$i=Object.assign({},Ki,{firstChunkSize:4e4,gps:[1,2,3,4]});async function is(n){let t=new Ue($i);await t.read(n);let e=await t.parse();if(e&&e.gps){let{latitude:i,longitude:r}=e.gps;return{latitude:i,longitude:r}}}const Ji=Object.assign({},Ki,{tiff:!1,ifd1:!0,mergeOutput:!1});async function rs(n){let t=new Ue(Ji);await t.read(n);let e=await t.extractThumbnail();return e&&ai?si.from(e):e}async function ss(n){let t=await this.thumbnail(n);if(t!==void 0){let e=new Blob([t]);return URL.createObjectURL(e)}}const Zi=Object.assign({},Ki,{firstChunkSize:4e4,ifd0:[274]});async function tr(n){let t=new Ue(Zi);await t.read(n);let e=await t.parse();if(e&&e.ifd0)return e.ifd0[274]}const er=Object.freeze({1:{dimensionSwapped:!1,scaleX:1,scaleY:1,deg:0,rad:0},2:{dimensionSwapped:!1,scaleX:-1,scaleY:1,deg:0,rad:0},3:{dimensionSwapped:!1,scaleX:1,scaleY:1,deg:180,rad:180*Math.PI/180},4:{dimensionSwapped:!1,scaleX:-1,scaleY:1,deg:180,rad:180*Math.PI/180},5:{dimensionSwapped:!0,scaleX:1,scaleY:-1,deg:90,rad:90*Math.PI/180},6:{dimensionSwapped:!0,scaleX:1,scaleY:1,deg:90,rad:90*Math.PI/180},7:{dimensionSwapped:!0,scaleX:1,scaleY:-1,deg:270,rad:270*Math.PI/180},8:{dimensionSwapped:!0,scaleX:1,scaleY:1,deg:270,rad:270*Math.PI/180}});let $e=!0,Je=!0;if(typeof navigator=="object"){let n=navigator.userAgent;if(n.includes("iPad")||n.includes("iPhone")){let t=n.match(/OS (\d+)_(\d+)/);if(t){let[,e,i]=t;$e=Number(e)+.1*Number(i)<13.4,Je=!1}}else if(n.includes("OS X 10")){let[,t]=n.match(/OS X 10[_.](\d+)/);$e=Je=Number(t)<15}if(n.includes("Chrome/")){let[,t]=n.match(/Chrome\/(\d+)/);$e=Je=Number(t)<81}else if(n.includes("Firefox/")){let[,t]=n.match(/Firefox\/(\d+)/);$e=Je=Number(t)<77}}async function as(n){let t=await tr(n);return Object.assign({canvas:$e,css:Je},er[t])}class Aa extends Wt{constructor(...t){super(...t),st(this,"ranges",new Ta),this.byteLength!==0&&this.ranges.add(0,this.byteLength)}_tryExtend(t,e,i){if(t===0&&this.byteLength===0&&i){let r=new DataView(i.buffer||i,i.byteOffset,i.byteLength);this._swapDataView(r)}else{let r=t+e;if(r>this.byteLength){let{dataView:s}=this._extend(r);this._swapDataView(s)}}}_extend(t){let e;e=ai?si.allocUnsafe(t):new Uint8Array(t);let i=new DataView(e.buffer,e.byteOffset,e.byteLength);return e.set(new Uint8Array(this.buffer,this.byteOffset,this.byteLength),0),{uintView:e,dataView:i}}subarray(t,e,i=!1){return e=e||this._lengthToEnd(t),i&&this._tryExtend(t,e),this.ranges.add(t,e),super.subarray(t,e)}set(t,e,i=!1){i&&this._tryExtend(e,t.byteLength,t);let r=super.set(t,e);return this.ranges.add(e,r.byteLength),r}async ensureChunk(t,e){this.chunked&&(this.ranges.available(t,e)||await this.readChunk(t,e))}available(t,e){return this.ranges.available(t,e)}}class Ta{constructor(){st(this,"list",[])}get length(){return this.list.length}add(t,e,i=0){let r=t+e,s=this.list.filter(o=>kr(t,o.offset,r)||kr(t,o.end,r));if(s.length>0){t=Math.min(t,...s.map(a=>a.offset)),r=Math.max(r,...s.map(a=>a.end)),e=r-t;let o=s.shift();o.offset=t,o.length=e,o.end=r,this.list=this.list.filter(a=>!s.includes(a))}else this.list.push({offset:t,length:e,end:r})}available(t,e){let i=t+e;return this.list.some(r=>r.offset<=t&&i<=r.end)}}function kr(n,t,e){return n<=t&&t<=e}class ci extends Aa{constructor(t,e){super(0),st(this,"chunksRead",0),this.input=t,this.options=e}async readWhole(){this.chunked=!1,await this.readChunk(this.nextChunkOffset)}async readChunked(){this.chunked=!0,await this.readChunk(0,this.options.firstChunkSize)}async readNextChunk(t=this.nextChunkOffset){if(this.fullyRead)return this.chunksRead++,!1;let e=this.options.chunkSize,i=await this.readChunk(t,e);return!!i&&i.byteLength===e}async readChunk(t,e){if(this.chunksRead++,(e=this.safeWrapAddress(t,e))!==0)return this._readChunk(t,e)}safeWrapAddress(t,e){return this.size!==void 0&&t+e>this.size?Math.max(0,this.size-t):e}get nextChunkOffset(){if(this.ranges.list.length!==0)return this.ranges.list[0].length}get canReadNextChunk(){return this.chunksRead<this.options.chunkLimit}get fullyRead(){return this.size!==void 0&&this.nextChunkOffset===this.size}read(){return this.options.chunked?this.readChunked():this.readWhole()}close(){}}be.set("blob",class extends ci{async readWhole(){this.chunked=!1;let n=await sn(this.input);this._swapArrayBuffer(n)}readChunked(){return this.chunked=!0,this.size=this.input.size,super.readChunked()}async _readChunk(n,t){let e=t?n+t:void 0,i=this.input.slice(n,e),r=await sn(i);return this.set(r,n,!0)}});var Ia=Object.freeze({__proto__:null,default:Fa,Exifr:Ue,fileParsers:ye,segmentParsers:Pt,fileReaders:be,tagKeys:Nt,tagValues:oe,tagRevivers:Ee,createDictionary:Ct,extendDictionary:an,fetchUrlAsArrayBuffer:rn,readBlobAsArrayBuffer:sn,chunkedProps:Te,otherSegments:gn,segments:on,tiffBlocks:St,segmentsAndBlocks:Ie,tiffExtractables:Ne,inheritables:yn,allFormatters:Re,Options:ln,parse:oi,gpsOnlyOptions:$i,gps:is,thumbnailOnlyOptions:Ji,thumbnail:rs,thumbnailUrl:ss,orientationOnlyOptions:Zi,orientation:tr,rotations:er,get rotateCanvas(){return $e},get rotateCss(){return Je},rotation:as});be.set("url",class extends ci{async readWhole(){this.chunked=!1;let n=await rn(this.input);n instanceof ArrayBuffer?this._swapArrayBuffer(n):n instanceof Uint8Array&&this._swapBuffer(n)}async _readChunk(n,t){let e=t?n+t-1:void 0,i=this.options.httpHeaders||{};(n||e)&&(i.range=`bytes=${[n,e].join("-")}`);let r=await Qi(this.input,{headers:i}),s=await r.arrayBuffer(),o=s.byteLength;if(r.status!==416)return o!==t&&(this.size=n+o),this.set(s,n,!0)}});Wt.prototype.getUint64=function(n){let t=this.getUint32(n),e=this.getUint32(n+4);return t<1048575?t<<32|e:typeof Vn!==void 0?(console.warn("Using BigInt because of type 64uint but JS can only handle 53b numbers."),Vn(t)<<Vn(32)|Vn(e)):void kt("Trying to read 64b value but JS can only handle 53b numbers.")};class Na extends li{parseBoxes(t=0){let e=[];for(;t<this.file.byteLength-4;){let i=this.parseBoxHead(t);if(e.push(i),i.length===0)break;t+=i.length}return e}parseSubBoxes(t){t.boxes=this.parseBoxes(t.start)}findBox(t,e){return t.boxes===void 0&&this.parseSubBoxes(t),t.boxes.find(i=>i.kind===e)}parseBoxHead(t){let e=this.file.getUint32(t),i=this.file.getString(t+4,4),r=t+8;return e===1&&(e=this.file.getUint64(t+8),r+=8),{offset:t,length:e,kind:i,start:r}}parseBoxFullHead(t){if(t.version!==void 0)return;let e=this.file.getUint32(t.start);t.version=e>>24,t.start+=4}}class os extends Na{static canHandle(t,e){if(e!==0)return!1;let i=t.getUint16(2);if(i>50)return!1;let r=16,s=[];for(;r<i;)s.push(t.getString(r,4)),r+=4;return s.includes(this.type)}async parse(){let t=this.file.getUint32(0),e=this.parseBoxHead(t);for(;e.kind!=="meta";)t+=e.length,await this.file.ensureChunk(t,16),e=this.parseBoxHead(t);await this.file.ensureChunk(e.offset,e.length),this.parseBoxFullHead(e),this.parseSubBoxes(e),this.options.icc.enabled&&await this.findIcc(e),this.options.tiff.enabled&&await this.findExif(e)}async registerSegment(t,e,i){await this.file.ensureChunk(e,i);let r=this.file.subarray(e,i);this.createParser(t,r)}async findIcc(t){let e=this.findBox(t,"iprp");if(e===void 0)return;let i=this.findBox(e,"ipco");if(i===void 0)return;let r=this.findBox(i,"colr");r!==void 0&&await this.registerSegment("icc",r.offset+12,r.length)}async findExif(t){let e=this.findBox(t,"iinf");if(e===void 0)return;let i=this.findBox(t,"iloc");if(i===void 0)return;let r=this.findExifLocIdInIinf(e),s=this.findExtentInIloc(i,r);if(s===void 0)return;let[o,a]=s;await this.file.ensureChunk(o,a);let l=4+this.file.getUint32(o);o+=l,a-=l,await this.registerSegment("tiff",o,a)}findExifLocIdInIinf(t){this.parseBoxFullHead(t);let e,i,r,s,o=t.start,a=this.file.getUint16(o);for(o+=2;a--;){if(e=this.parseBoxHead(o),this.parseBoxFullHead(e),i=e.start,e.version>=2&&(r=e.version===3?4:2,s=this.file.getString(i+r+2,4),s==="Exif"))return this.file.getUintBytes(i,r);o+=e.length}}get8bits(t){let e=this.file.getUint8(t);return[e>>4,15&e]}findExtentInIloc(t,e){this.parseBoxFullHead(t);let i=t.start,[r,s]=this.get8bits(i++),[o,a]=this.get8bits(i++),l=t.version===2?4:2,u=t.version===1||t.version===2?2:0,h=a+r+s,c=t.version===2?4:2,m=this.file.getUintBytes(i,c);for(i+=c;m--;){let d=this.file.getUintBytes(i,l);i+=l+u+2+o;let f=this.file.getUint16(i);if(i+=2,d===e)return f>1&&console.warn(`ILOC box has more than one extent but we're only processing one
Please create an issue at https://github.com/MikeKovarik/exifr with this file`),[this.file.getUintBytes(i+a,r),this.file.getUintBytes(i+a+r,s)];i+=f*h}}}class ls extends os{}st(ls,"type","heic");class Fr extends os{}st(Fr,"type","avif"),ye.set("heic",ls),ye.set("avif",Fr),Ct(Nt,["ifd0","ifd1"],[[256,"ImageWidth"],[257,"ImageHeight"],[258,"BitsPerSample"],[259,"Compression"],[262,"PhotometricInterpretation"],[270,"ImageDescription"],[271,"Make"],[272,"Model"],[273,"StripOffsets"],[274,"Orientation"],[277,"SamplesPerPixel"],[278,"RowsPerStrip"],[279,"StripByteCounts"],[282,"XResolution"],[283,"YResolution"],[284,"PlanarConfiguration"],[296,"ResolutionUnit"],[301,"TransferFunction"],[305,"Software"],[306,"ModifyDate"],[315,"Artist"],[316,"HostComputer"],[317,"Predictor"],[318,"WhitePoint"],[319,"PrimaryChromaticities"],[513,"ThumbnailOffset"],[514,"ThumbnailLength"],[529,"YCbCrCoefficients"],[530,"YCbCrSubSampling"],[531,"YCbCrPositioning"],[532,"ReferenceBlackWhite"],[700,"ApplicationNotes"],[33432,"Copyright"],[33723,"IPTC"],[34665,"ExifIFD"],[34675,"ICC"],[34853,"GpsIFD"],[330,"SubIFD"],[40965,"InteropIFD"],[40091,"XPTitle"],[40092,"XPComment"],[40093,"XPAuthor"],[40094,"XPKeywords"],[40095,"XPSubject"]]),Ct(Nt,"exif",[[33434,"ExposureTime"],[33437,"FNumber"],[34850,"ExposureProgram"],[34852,"SpectralSensitivity"],[34855,"ISO"],[34858,"TimeZoneOffset"],[34859,"SelfTimerMode"],[34864,"SensitivityType"],[34865,"StandardOutputSensitivity"],[34866,"RecommendedExposureIndex"],[34867,"ISOSpeed"],[34868,"ISOSpeedLatitudeyyy"],[34869,"ISOSpeedLatitudezzz"],[36864,"ExifVersion"],[36867,"DateTimeOriginal"],[36868,"CreateDate"],[36873,"GooglePlusUploadCode"],[36880,"OffsetTime"],[36881,"OffsetTimeOriginal"],[36882,"OffsetTimeDigitized"],[37121,"ComponentsConfiguration"],[37122,"CompressedBitsPerPixel"],[37377,"ShutterSpeedValue"],[37378,"ApertureValue"],[37379,"BrightnessValue"],[37380,"ExposureCompensation"],[37381,"MaxApertureValue"],[37382,"SubjectDistance"],[37383,"MeteringMode"],[37384,"LightSource"],[37385,"Flash"],[37386,"FocalLength"],[37393,"ImageNumber"],[37394,"SecurityClassification"],[37395,"ImageHistory"],[37396,"SubjectArea"],[37500,"MakerNote"],[37510,"UserComment"],[37520,"SubSecTime"],[37521,"SubSecTimeOriginal"],[37522,"SubSecTimeDigitized"],[37888,"AmbientTemperature"],[37889,"Humidity"],[37890,"Pressure"],[37891,"WaterDepth"],[37892,"Acceleration"],[37893,"CameraElevationAngle"],[40960,"FlashpixVersion"],[40961,"ColorSpace"],[40962,"ExifImageWidth"],[40963,"ExifImageHeight"],[40964,"RelatedSoundFile"],[41483,"FlashEnergy"],[41486,"FocalPlaneXResolution"],[41487,"FocalPlaneYResolution"],[41488,"FocalPlaneResolutionUnit"],[41492,"SubjectLocation"],[41493,"ExposureIndex"],[41495,"SensingMethod"],[41728,"FileSource"],[41729,"SceneType"],[41730,"CFAPattern"],[41985,"CustomRendered"],[41986,"ExposureMode"],[41987,"WhiteBalance"],[41988,"DigitalZoomRatio"],[41989,"FocalLengthIn35mmFormat"],[41990,"SceneCaptureType"],[41991,"GainControl"],[41992,"Contrast"],[41993,"Saturation"],[41994,"Sharpness"],[41996,"SubjectDistanceRange"],[42016,"ImageUniqueID"],[42032,"OwnerName"],[42033,"SerialNumber"],[42034,"LensInfo"],[42035,"LensMake"],[42036,"LensModel"],[42037,"LensSerialNumber"],[42080,"CompositeImage"],[42081,"CompositeImageCount"],[42082,"CompositeImageExposureTimes"],[42240,"Gamma"],[59932,"Padding"],[59933,"OffsetSchema"],[65e3,"OwnerName"],[65001,"SerialNumber"],[65002,"Lens"],[65100,"RawFile"],[65101,"Converter"],[65102,"WhiteBalance"],[65105,"Exposure"],[65106,"Shadows"],[65107,"Brightness"],[65108,"Contrast"],[65109,"Saturation"],[65110,"Sharpness"],[65111,"Smoothness"],[65112,"MoireFilter"],[40965,"InteropIFD"]]),Ct(Nt,"gps",[[0,"GPSVersionID"],[1,"GPSLatitudeRef"],[2,"GPSLatitude"],[3,"GPSLongitudeRef"],[4,"GPSLongitude"],[5,"GPSAltitudeRef"],[6,"GPSAltitude"],[7,"GPSTimeStamp"],[8,"GPSSatellites"],[9,"GPSStatus"],[10,"GPSMeasureMode"],[11,"GPSDOP"],[12,"GPSSpeedRef"],[13,"GPSSpeed"],[14,"GPSTrackRef"],[15,"GPSTrack"],[16,"GPSImgDirectionRef"],[17,"GPSImgDirection"],[18,"GPSMapDatum"],[19,"GPSDestLatitudeRef"],[20,"GPSDestLatitude"],[21,"GPSDestLongitudeRef"],[22,"GPSDestLongitude"],[23,"GPSDestBearingRef"],[24,"GPSDestBearing"],[25,"GPSDestDistanceRef"],[26,"GPSDestDistance"],[27,"GPSProcessingMethod"],[28,"GPSAreaInformation"],[29,"GPSDateStamp"],[30,"GPSDifferential"],[31,"GPSHPositioningError"]]),Ct(oe,["ifd0","ifd1"],[[274,{1:"Horizontal (normal)",2:"Mirror horizontal",3:"Rotate 180",4:"Mirror vertical",5:"Mirror horizontal and rotate 270 CW",6:"Rotate 90 CW",7:"Mirror horizontal and rotate 90 CW",8:"Rotate 270 CW"}],[296,{1:"None",2:"inches",3:"cm"}]]);let Pn=Ct(oe,"exif",[[34850,{0:"Not defined",1:"Manual",2:"Normal program",3:"Aperture priority",4:"Shutter priority",5:"Creative program",6:"Action program",7:"Portrait mode",8:"Landscape mode"}],[37121,{0:"-",1:"Y",2:"Cb",3:"Cr",4:"R",5:"G",6:"B"}],[37383,{0:"Unknown",1:"Average",2:"CenterWeightedAverage",3:"Spot",4:"MultiSpot",5:"Pattern",6:"Partial",255:"Other"}],[37384,{0:"Unknown",1:"Daylight",2:"Fluorescent",3:"Tungsten (incandescent light)",4:"Flash",9:"Fine weather",10:"Cloudy weather",11:"Shade",12:"Daylight fluorescent (D 5700 - 7100K)",13:"Day white fluorescent (N 4600 - 5400K)",14:"Cool white fluorescent (W 3900 - 4500K)",15:"White fluorescent (WW 3200 - 3700K)",17:"Standard light A",18:"Standard light B",19:"Standard light C",20:"D55",21:"D65",22:"D75",23:"D50",24:"ISO studio tungsten",255:"Other"}],[37385,{0:"Flash did not fire",1:"Flash fired",5:"Strobe return light not detected",7:"Strobe return light detected",9:"Flash fired, compulsory flash mode",13:"Flash fired, compulsory flash mode, return light not detected",15:"Flash fired, compulsory flash mode, return light detected",16:"Flash did not fire, compulsory flash mode",24:"Flash did not fire, auto mode",25:"Flash fired, auto mode",29:"Flash fired, auto mode, return light not detected",31:"Flash fired, auto mode, return light detected",32:"No flash function",65:"Flash fired, red-eye reduction mode",69:"Flash fired, red-eye reduction mode, return light not detected",71:"Flash fired, red-eye reduction mode, return light detected",73:"Flash fired, compulsory flash mode, red-eye reduction mode",77:"Flash fired, compulsory flash mode, red-eye reduction mode, return light not detected",79:"Flash fired, compulsory flash mode, red-eye reduction mode, return light detected",89:"Flash fired, auto mode, red-eye reduction mode",93:"Flash fired, auto mode, return light not detected, red-eye reduction mode",95:"Flash fired, auto mode, return light detected, red-eye reduction mode"}],[41495,{1:"Not defined",2:"One-chip color area sensor",3:"Two-chip color area sensor",4:"Three-chip color area sensor",5:"Color sequential area sensor",7:"Trilinear sensor",8:"Color sequential linear sensor"}],[41728,{1:"Film Scanner",2:"Reflection Print Scanner",3:"Digital Camera"}],[41729,{1:"Directly photographed"}],[41985,{0:"Normal",1:"Custom",2:"HDR (no original saved)",3:"HDR (original saved)",4:"Original (for HDR)",6:"Panorama",7:"Portrait HDR",8:"Portrait"}],[41986,{0:"Auto",1:"Manual",2:"Auto bracket"}],[41987,{0:"Auto",1:"Manual"}],[41990,{0:"Standard",1:"Landscape",2:"Portrait",3:"Night",4:"Other"}],[41991,{0:"None",1:"Low gain up",2:"High gain up",3:"Low gain down",4:"High gain down"}],[41996,{0:"Unknown",1:"Macro",2:"Close",3:"Distant"}],[42080,{0:"Unknown",1:"Not a Composite Image",2:"General Composite Image",3:"Composite Image Captured While Shooting"}]]);const Ar={1:"No absolute unit of measurement",2:"Inch",3:"Centimeter"};Pn.set(37392,Ar),Pn.set(41488,Ar);const Ni={0:"Normal",1:"Low",2:"High"};function Tr(n){return typeof n=="object"&&n.length!==void 0?n[0]:n}function Ir(n){let t=Array.from(n).slice(1);return t[1]>15&&(t=t.map(e=>String.fromCharCode(e))),t[2]!=="0"&&t[2]!==0||t.pop(),t.join(".")}function Ri(n){if(typeof n=="string"){var[t,e,i,r,s,o]=n.trim().split(/[-: ]/g).map(Number),a=new Date(t,e-1,i);return Number.isNaN(r)||Number.isNaN(s)||Number.isNaN(o)||(a.setHours(r),a.setMinutes(s),a.setSeconds(o)),Number.isNaN(+a)?n:a}}function vn(n){if(typeof n=="string")return n;let t=[];if(n[1]===0&&n[n.length-1]===0)for(let e=0;e<n.length;e+=2)t.push(Nr(n[e+1],n[e]));else for(let e=0;e<n.length;e+=2)t.push(Nr(n[e],n[e+1]));return Ke(String.fromCodePoint(...t))}function Nr(n,t){return n<<8|t}Pn.set(41992,Ni),Pn.set(41993,Ni),Pn.set(41994,Ni),Ct(Ee,["ifd0","ifd1"],[[50827,function(n){return typeof n!="string"?es(n):n}],[306,Ri],[40091,vn],[40092,vn],[40093,vn],[40094,vn],[40095,vn]]),Ct(Ee,"exif",[[40960,Ir],[36864,Ir],[36867,Ri],[36868,Ri],[40962,Tr],[40963,Tr]]),Ct(Ee,"gps",[[0,n=>Array.from(n).join(".")],[7,n=>Array.from(n).join(":")]]);class Li extends ae{static canHandle(t,e){return t.getUint8(e+1)===225&&t.getUint32(e+4)===1752462448&&t.getString(e+4,20)==="http://ns.adobe.com/"}static headerLength(t,e){return t.getString(e+4,34)==="http://ns.adobe.com/xmp/extension/"?79:33}static findPosition(t,e){let i=super.findPosition(t,e);return i.multiSegment=i.extended=i.headerLength===79,i.multiSegment?(i.chunkCount=t.getUint8(e+72),i.chunkNumber=t.getUint8(e+76),t.getUint8(e+77)!==0&&i.chunkNumber++):(i.chunkCount=1/0,i.chunkNumber=-1),i}static handleMultiSegments(t){return t.map(e=>e.chunk.getString()).join("")}normalizeInput(t){return typeof t=="string"?t:Wt.from(t).getString()}parse(t=this.chunk){if(!this.localOptions.parse)return t;t=function(s){let o={},a={};for(let l of ds)o[l]=[],a[l]=0;return s.replace(Ua,(l,u,h)=>{if(u==="<"){let c=++a[h];return o[h].push(c),`${l}#${c}`}return`${l}#${o[h].pop()}`})}(t);let e=dn.findAll(t,"rdf","Description");e.length===0&&e.push(new dn("rdf","Description",void 0,t));let i,r={};for(let s of e)for(let o of s.properties)i=Ea(o.ns,r),cs(o,i);return function(s){let o;for(let a in s)o=s[a]=Jn(s[a]),o===void 0&&delete s[a];return Jn(s)}(r)}assignToOutput(t,e){if(this.localOptions.parse)for(let[i,r]of Object.entries(e))switch(i){case"tiff":this.assignObjectToOutput(t,"ifd0",r);break;case"exif":this.assignObjectToOutput(t,"exif",r);break;case"xmlns":break;default:this.assignObjectToOutput(t,i,r)}else t.xmp=e}}st(Li,"type","xmp"),st(Li,"multiSegment",!0),Pt.set("xmp",Li);class ti{static findAll(t){return us(t,/([a-zA-Z0-9-]+):([a-zA-Z0-9-]+)=("[^"]*"|'[^']*')/gm).map(ti.unpackMatch)}static unpackMatch(t){let e=t[1],i=t[2],r=t[3].slice(1,-1);return r=hs(r),new ti(e,i,r)}constructor(t,e,i){this.ns=t,this.name=e,this.value=i}serialize(){return this.value}}class dn{static findAll(t,e,i){if(e!==void 0||i!==void 0){e=e||"[\\w\\d-]+",i=i||"[\\w\\d-]+";var r=new RegExp(`<(${e}):(${i})(#\\d+)?((\\s+?[\\w\\d-:]+=("[^"]*"|'[^']*'))*\\s*)(\\/>|>([\\s\\S]*?)<\\/\\1:\\2\\3>)`,"gm")}else r=/<([\w\d-]+):([\w\d-]+)(#\d+)?((\s+?[\w\d-:]+=("[^"]*"|'[^']*'))*\s*)(\/>|>([\s\S]*?)<\/\1:\2\3>)/gm;return us(t,r).map(dn.unpackMatch)}static unpackMatch(t){let e=t[1],i=t[2],r=t[4],s=t[8];return new dn(e,i,r,s)}constructor(t,e,i,r){this.ns=t,this.name=e,this.attrString=i,this.innerXml=r,this.attrs=ti.findAll(i),this.children=dn.findAll(r),this.value=this.children.length===0?hs(r):void 0,this.properties=[...this.attrs,...this.children]}get isPrimitive(){return this.value!==void 0&&this.attrs.length===0&&this.children.length===0}get isListContainer(){return this.children.length===1&&this.children[0].isList}get isList(){let{ns:t,name:e}=this;return t==="rdf"&&(e==="Seq"||e==="Bag"||e==="Alt")}get isListItem(){return this.ns==="rdf"&&this.name==="li"}serialize(){if(this.properties.length===0&&this.value===void 0)return;if(this.isPrimitive)return this.value;if(this.isListContainer)return this.children[0].serialize();if(this.isList)return La(this.children.map(Ra));if(this.isListItem&&this.children.length===1&&this.attrs.length===0)return this.children[0].serialize();let t={};for(let e of this.properties)cs(e,t);return this.value!==void 0&&(t.value=this.value),Jn(t)}}function cs(n,t){let e=n.serialize();e!==void 0&&(t[n.name]=e)}var Ra=n=>n.serialize(),La=n=>n.length===1?n[0]:n,Ea=(n,t)=>t[n]?t[n]:t[n]={};function us(n,t){let e,i=[];if(!n)return i;for(;(e=t.exec(n))!==null;)i.push(e);return i}function hs(n){if(function(i){return i==null||i==="null"||i==="undefined"||i===""||i.trim()===""}(n))return;let t=Number(n);if(!Number.isNaN(t))return t;let e=n.toLowerCase();return e==="true"||e!=="false"&&n.trim()}const ds=["rdf:li","rdf:Seq","rdf:Bag","rdf:Alt","rdf:Description"],Ua=new RegExp(`(<|\\/)(${ds.join("|")})`,"g");var Ba=Object.freeze({__proto__:null,default:Ia,Exifr:Ue,fileParsers:ye,segmentParsers:Pt,fileReaders:be,tagKeys:Nt,tagValues:oe,tagRevivers:Ee,createDictionary:Ct,extendDictionary:an,fetchUrlAsArrayBuffer:rn,readBlobAsArrayBuffer:sn,chunkedProps:Te,otherSegments:gn,segments:on,tiffBlocks:St,segmentsAndBlocks:Ie,tiffExtractables:Ne,inheritables:yn,allFormatters:Re,Options:ln,parse:oi,gpsOnlyOptions:$i,gps:is,thumbnailOnlyOptions:Ji,thumbnail:rs,thumbnailUrl:ss,orientationOnlyOptions:Zi,orientation:tr,rotations:er,get rotateCanvas(){return $e},get rotateCss(){return Je},rotation:as});let Rr=$n("fs",n=>n.promises);be.set("fs",class extends ci{async readWhole(){this.chunked=!1,this.fs=await Rr;let n=await this.fs.readFile(this.input);this._swapBuffer(n)}async readChunked(){this.chunked=!0,this.fs=await Rr,await this.open(),await this.readChunk(0,this.options.firstChunkSize)}async open(){this.fh===void 0&&(this.fh=await this.fs.open(this.input,"r"),this.size=(await this.fh.stat(this.input)).size)}async _readChunk(n,t){this.fh===void 0&&await this.open(),n+t>this.size&&(t=this.size-n);var e=this.subarray(n,t,!0);return await this.fh.read(e.dataView,0,t,n),e}async close(){if(this.fh){let n=this.fh;this.fh=void 0,await n.close()}}});be.set("base64",class extends ci{constructor(...n){super(...n),this.input=this.input.replace(/^data:([^;]+);base64,/gim,""),this.size=this.input.length/4*3,this.input.endsWith("==")?this.size-=2:this.input.endsWith("=")&&(this.size-=1)}async _readChunk(n,t){let e,i,r=this.input;n===void 0?(n=0,e=0,i=0):(e=4*Math.floor(n/3),i=n-e/4*3),t===void 0&&(t=this.size);let s=n+t,o=e+4*Math.ceil(s/3);r=r.slice(e,o);let a=Math.min(t,this.size-n);if(ai){let l=si.from(r,"base64").slice(i,i+a);return this.set(l,n,!0)}{let l=this.subarray(n,a,!0),u=atob(r),h=l.toUint8();for(let c=0;c<a;c++)h[c]=u.charCodeAt(i+c);return l}}});class Lr extends li{static canHandle(t,e){return e===18761||e===19789}extendOptions(t){let{ifd0:e,xmp:i,iptc:r,icc:s}=t;i.enabled&&e.deps.add(700),r.enabled&&e.deps.add(33723),s.enabled&&e.deps.add(34675),e.finalizeFilters()}async parse(){let{tiff:t,xmp:e,iptc:i,icc:r}=this.options;if(t.enabled||e.enabled||i.enabled||r.enabled){let s=Math.max(Gi(this.options),this.options.chunkSize);await this.file.ensureChunk(0,s),this.createParser("tiff",this.file),this.parsers.tiff.parseHeader(),await this.parsers.tiff.parseIfd0Block(),this.adaptTiffPropAsSegment("xmp"),this.adaptTiffPropAsSegment("iptc"),this.adaptTiffPropAsSegment("icc")}}adaptTiffPropAsSegment(t){if(this.parsers.tiff[t]){let e=this.parsers.tiff[t];this.injectSegment(t,e)}}}st(Lr,"type","tiff"),ye.set("tiff",Lr);let Da=$n("zlib");const Oa=["ihdr","iccp","text","itxt","exif"];class Er extends li{constructor(...t){super(...t),st(this,"catchError",e=>this.errors.push(e)),st(this,"metaChunks",[]),st(this,"unknownChunks",[])}static canHandle(t,e){return e===35152&&t.getUint32(0)===2303741511&&t.getUint32(4)===218765834}async parse(){let{file:t}=this;await this.findPngChunksInRange(8,t.byteLength),await this.readSegments(this.metaChunks),this.findIhdr(),this.parseTextChunks(),await this.findExif().catch(this.catchError),await this.findXmp().catch(this.catchError),await this.findIcc().catch(this.catchError)}async findPngChunksInRange(t,e){let{file:i}=this;for(;t<e;){let r=i.getUint32(t),s=i.getUint32(t+4),o=i.getString(t+4,4).toLowerCase(),a=r+4+4+4,l={type:o,offset:t,length:a,start:t+4+4,size:r,marker:s};Oa.includes(o)?this.metaChunks.push(l):this.unknownChunks.push(l),t+=a}}parseTextChunks(){let t=this.metaChunks.filter(e=>e.type==="text");for(let e of t){let[i,r]=this.file.getString(e.start,e.size).split("\0");this.injectKeyValToIhdr(i,r)}}injectKeyValToIhdr(t,e){let i=this.parsers.ihdr;i&&i.raw.set(t,e)}findIhdr(){let t=this.metaChunks.find(e=>e.type==="ihdr");t&&this.options.ihdr.enabled!==!1&&this.createParser("ihdr",t.chunk)}async findExif(){let t=this.metaChunks.find(e=>e.type==="exif");t&&this.injectSegment("tiff",t.chunk)}async findXmp(){let t=this.metaChunks.filter(e=>e.type==="itxt");for(let e of t)e.chunk.getString(0,17)==="XML:com.adobe.xmp"&&this.injectSegment("xmp",e.chunk)}async findIcc(){let t=this.metaChunks.find(a=>a.type==="iccp");if(!t)return;let{chunk:e}=t,i=e.getUint8Array(0,81),r=0;for(;r<80&&i[r]!==0;)r++;let s=r+2,o=e.getString(0,r);if(this.injectKeyValToIhdr("ProfileName",o),Kn){let a=await Da,l=e.getUint8Array(s);l=a.inflateSync(l),this.injectSegment("icc",l)}}}st(Er,"type","png"),ye.set("png",Er),Ct(Nt,"interop",[[1,"InteropIndex"],[2,"InteropVersion"],[4096,"RelatedImageFileFormat"],[4097,"RelatedImageWidth"],[4098,"RelatedImageHeight"]]),an(Nt,"ifd0",[[11,"ProcessingSoftware"],[254,"SubfileType"],[255,"OldSubfileType"],[263,"Thresholding"],[264,"CellWidth"],[265,"CellLength"],[266,"FillOrder"],[269,"DocumentName"],[280,"MinSampleValue"],[281,"MaxSampleValue"],[285,"PageName"],[286,"XPosition"],[287,"YPosition"],[290,"GrayResponseUnit"],[297,"PageNumber"],[321,"HalftoneHints"],[322,"TileWidth"],[323,"TileLength"],[332,"InkSet"],[337,"TargetPrinter"],[18246,"Rating"],[18249,"RatingPercent"],[33550,"PixelScale"],[34264,"ModelTransform"],[34377,"PhotoshopSettings"],[50706,"DNGVersion"],[50707,"DNGBackwardVersion"],[50708,"UniqueCameraModel"],[50709,"LocalizedCameraModel"],[50736,"DNGLensInfo"],[50739,"ShadowScale"],[50740,"DNGPrivateData"],[33920,"IntergraphMatrix"],[33922,"ModelTiePoint"],[34118,"SEMInfo"],[34735,"GeoTiffDirectory"],[34736,"GeoTiffDoubleParams"],[34737,"GeoTiffAsciiParams"],[50341,"PrintIM"],[50721,"ColorMatrix1"],[50722,"ColorMatrix2"],[50723,"CameraCalibration1"],[50724,"CameraCalibration2"],[50725,"ReductionMatrix1"],[50726,"ReductionMatrix2"],[50727,"AnalogBalance"],[50728,"AsShotNeutral"],[50729,"AsShotWhiteXY"],[50730,"BaselineExposure"],[50731,"BaselineNoise"],[50732,"BaselineSharpness"],[50734,"LinearResponseLimit"],[50735,"CameraSerialNumber"],[50741,"MakerNoteSafety"],[50778,"CalibrationIlluminant1"],[50779,"CalibrationIlluminant2"],[50781,"RawDataUniqueID"],[50827,"OriginalRawFileName"],[50828,"OriginalRawFileData"],[50831,"AsShotICCProfile"],[50832,"AsShotPreProfileMatrix"],[50833,"CurrentICCProfile"],[50834,"CurrentPreProfileMatrix"],[50879,"ColorimetricReference"],[50885,"SRawType"],[50898,"PanasonicTitle"],[50899,"PanasonicTitle2"],[50931,"CameraCalibrationSig"],[50932,"ProfileCalibrationSig"],[50933,"ProfileIFD"],[50934,"AsShotProfileName"],[50936,"ProfileName"],[50937,"ProfileHueSatMapDims"],[50938,"ProfileHueSatMapData1"],[50939,"ProfileHueSatMapData2"],[50940,"ProfileToneCurve"],[50941,"ProfileEmbedPolicy"],[50942,"ProfileCopyright"],[50964,"ForwardMatrix1"],[50965,"ForwardMatrix2"],[50966,"PreviewApplicationName"],[50967,"PreviewApplicationVersion"],[50968,"PreviewSettingsName"],[50969,"PreviewSettingsDigest"],[50970,"PreviewColorSpace"],[50971,"PreviewDateTime"],[50972,"RawImageDigest"],[50973,"OriginalRawFileDigest"],[50981,"ProfileLookTableDims"],[50982,"ProfileLookTableData"],[51043,"TimeCodes"],[51044,"FrameRate"],[51058,"TStop"],[51081,"ReelName"],[51089,"OriginalDefaultFinalSize"],[51090,"OriginalBestQualitySize"],[51091,"OriginalDefaultCropSize"],[51105,"CameraLabel"],[51107,"ProfileHueSatMapEncoding"],[51108,"ProfileLookTableEncoding"],[51109,"BaselineExposureOffset"],[51110,"DefaultBlackRender"],[51111,"NewRawImageDigest"],[51112,"RawToPreviewGain"]]);let Ur=[[273,"StripOffsets"],[279,"StripByteCounts"],[288,"FreeOffsets"],[289,"FreeByteCounts"],[291,"GrayResponseCurve"],[292,"T4Options"],[293,"T6Options"],[300,"ColorResponseUnit"],[320,"ColorMap"],[324,"TileOffsets"],[325,"TileByteCounts"],[326,"BadFaxLines"],[327,"CleanFaxData"],[328,"ConsecutiveBadFaxLines"],[330,"SubIFD"],[333,"InkNames"],[334,"NumberofInks"],[336,"DotRange"],[338,"ExtraSamples"],[339,"SampleFormat"],[340,"SMinSampleValue"],[341,"SMaxSampleValue"],[342,"TransferRange"],[343,"ClipPath"],[344,"XClipPathUnits"],[345,"YClipPathUnits"],[346,"Indexed"],[347,"JPEGTables"],[351,"OPIProxy"],[400,"GlobalParametersIFD"],[401,"ProfileType"],[402,"FaxProfile"],[403,"CodingMethods"],[404,"VersionYear"],[405,"ModeNumber"],[433,"Decode"],[434,"DefaultImageColor"],[435,"T82Options"],[437,"JPEGTables"],[512,"JPEGProc"],[515,"JPEGRestartInterval"],[517,"JPEGLosslessPredictors"],[518,"JPEGPointTransforms"],[519,"JPEGQTables"],[520,"JPEGDCTables"],[521,"JPEGACTables"],[559,"StripRowCounts"],[999,"USPTOMiscellaneous"],[18247,"XP_DIP_XML"],[18248,"StitchInfo"],[28672,"SonyRawFileType"],[28688,"SonyToneCurve"],[28721,"VignettingCorrection"],[28722,"VignettingCorrParams"],[28724,"ChromaticAberrationCorrection"],[28725,"ChromaticAberrationCorrParams"],[28726,"DistortionCorrection"],[28727,"DistortionCorrParams"],[29895,"SonyCropTopLeft"],[29896,"SonyCropSize"],[32781,"ImageID"],[32931,"WangTag1"],[32932,"WangAnnotation"],[32933,"WangTag3"],[32934,"WangTag4"],[32953,"ImageReferencePoints"],[32954,"RegionXformTackPoint"],[32955,"WarpQuadrilateral"],[32956,"AffineTransformMat"],[32995,"Matteing"],[32996,"DataType"],[32997,"ImageDepth"],[32998,"TileDepth"],[33300,"ImageFullWidth"],[33301,"ImageFullHeight"],[33302,"TextureFormat"],[33303,"WrapModes"],[33304,"FovCot"],[33305,"MatrixWorldToScreen"],[33306,"MatrixWorldToCamera"],[33405,"Model2"],[33421,"CFARepeatPatternDim"],[33422,"CFAPattern2"],[33423,"BatteryLevel"],[33424,"KodakIFD"],[33445,"MDFileTag"],[33446,"MDScalePixel"],[33447,"MDColorTable"],[33448,"MDLabName"],[33449,"MDSampleInfo"],[33450,"MDPrepDate"],[33451,"MDPrepTime"],[33452,"MDFileUnits"],[33589,"AdventScale"],[33590,"AdventRevision"],[33628,"UIC1Tag"],[33629,"UIC2Tag"],[33630,"UIC3Tag"],[33631,"UIC4Tag"],[33918,"IntergraphPacketData"],[33919,"IntergraphFlagRegisters"],[33921,"INGRReserved"],[34016,"Site"],[34017,"ColorSequence"],[34018,"IT8Header"],[34019,"RasterPadding"],[34020,"BitsPerRunLength"],[34021,"BitsPerExtendedRunLength"],[34022,"ColorTable"],[34023,"ImageColorIndicator"],[34024,"BackgroundColorIndicator"],[34025,"ImageColorValue"],[34026,"BackgroundColorValue"],[34027,"PixelIntensityRange"],[34028,"TransparencyIndicator"],[34029,"ColorCharacterization"],[34030,"HCUsage"],[34031,"TrapIndicator"],[34032,"CMYKEquivalent"],[34152,"AFCP_IPTC"],[34232,"PixelMagicJBIGOptions"],[34263,"JPLCartoIFD"],[34306,"WB_GRGBLevels"],[34310,"LeafData"],[34687,"TIFF_FXExtensions"],[34688,"MultiProfiles"],[34689,"SharedData"],[34690,"T88Options"],[34732,"ImageLayer"],[34750,"JBIGOptions"],[34856,"Opto-ElectricConvFactor"],[34857,"Interlace"],[34908,"FaxRecvParams"],[34909,"FaxSubAddress"],[34910,"FaxRecvTime"],[34929,"FedexEDR"],[34954,"LeafSubIFD"],[37387,"FlashEnergy"],[37388,"SpatialFrequencyResponse"],[37389,"Noise"],[37390,"FocalPlaneXResolution"],[37391,"FocalPlaneYResolution"],[37392,"FocalPlaneResolutionUnit"],[37397,"ExposureIndex"],[37398,"TIFF-EPStandardID"],[37399,"SensingMethod"],[37434,"CIP3DataFile"],[37435,"CIP3Sheet"],[37436,"CIP3Side"],[37439,"StoNits"],[37679,"MSDocumentText"],[37680,"MSPropertySetStorage"],[37681,"MSDocumentTextPosition"],[37724,"ImageSourceData"],[40965,"InteropIFD"],[40976,"SamsungRawPointersOffset"],[40977,"SamsungRawPointersLength"],[41217,"SamsungRawByteOrder"],[41218,"SamsungRawUnknown"],[41484,"SpatialFrequencyResponse"],[41485,"Noise"],[41489,"ImageNumber"],[41490,"SecurityClassification"],[41491,"ImageHistory"],[41494,"TIFF-EPStandardID"],[41995,"DeviceSettingDescription"],[42112,"GDALMetadata"],[42113,"GDALNoData"],[44992,"ExpandSoftware"],[44993,"ExpandLens"],[44994,"ExpandFilm"],[44995,"ExpandFilterLens"],[44996,"ExpandScanner"],[44997,"ExpandFlashLamp"],[46275,"HasselbladRawImage"],[48129,"PixelFormat"],[48130,"Transformation"],[48131,"Uncompressed"],[48132,"ImageType"],[48256,"ImageWidth"],[48257,"ImageHeight"],[48258,"WidthResolution"],[48259,"HeightResolution"],[48320,"ImageOffset"],[48321,"ImageByteCount"],[48322,"AlphaOffset"],[48323,"AlphaByteCount"],[48324,"ImageDataDiscard"],[48325,"AlphaDataDiscard"],[50215,"OceScanjobDesc"],[50216,"OceApplicationSelector"],[50217,"OceIDNumber"],[50218,"OceImageLogic"],[50255,"Annotations"],[50459,"HasselbladExif"],[50547,"OriginalFileName"],[50560,"USPTOOriginalContentType"],[50656,"CR2CFAPattern"],[50710,"CFAPlaneColor"],[50711,"CFALayout"],[50712,"LinearizationTable"],[50713,"BlackLevelRepeatDim"],[50714,"BlackLevel"],[50715,"BlackLevelDeltaH"],[50716,"BlackLevelDeltaV"],[50717,"WhiteLevel"],[50718,"DefaultScale"],[50719,"DefaultCropOrigin"],[50720,"DefaultCropSize"],[50733,"BayerGreenSplit"],[50737,"ChromaBlurRadius"],[50738,"AntiAliasStrength"],[50752,"RawImageSegmentation"],[50780,"BestQualityScale"],[50784,"AliasLayerMetadata"],[50829,"ActiveArea"],[50830,"MaskedAreas"],[50935,"NoiseReductionApplied"],[50974,"SubTileBlockSize"],[50975,"RowInterleaveFactor"],[51008,"OpcodeList1"],[51009,"OpcodeList2"],[51022,"OpcodeList3"],[51041,"NoiseProfile"],[51114,"CacheVersion"],[51125,"DefaultUserCrop"],[51157,"NikonNEFInfo"],[65024,"KdcIFD"]];an(Nt,"ifd0",Ur),an(Nt,"exif",Ur),Ct(oe,"gps",[[23,{M:"Magnetic North",T:"True North"}],[25,{K:"Kilometers",M:"Miles",N:"Nautical Miles"}]]);class Ei extends ae{static canHandle(t,e){return t.getUint8(e+1)===224&&t.getUint32(e+4)===1246120262&&t.getUint8(e+8)===0}parse(){return this.parseTags(),this.translate(),this.output}parseTags(){this.raw=new Map([[0,this.chunk.getUint16(0)],[2,this.chunk.getUint8(2)],[3,this.chunk.getUint16(3)],[5,this.chunk.getUint16(5)],[7,this.chunk.getUint8(7)],[8,this.chunk.getUint8(8)]])}}st(Ei,"type","jfif"),st(Ei,"headerLength",9),Pt.set("jfif",Ei),Ct(Nt,"jfif",[[0,"JFIFVersion"],[2,"ResolutionUnit"],[3,"XResolution"],[5,"YResolution"],[7,"ThumbnailWidth"],[8,"ThumbnailHeight"]]);class Br extends ae{parse(){return this.parseTags(),this.translate(),this.output}parseTags(){this.raw=new Map([[0,this.chunk.getUint32(0)],[4,this.chunk.getUint32(4)],[8,this.chunk.getUint8(8)],[9,this.chunk.getUint8(9)],[10,this.chunk.getUint8(10)],[11,this.chunk.getUint8(11)],[12,this.chunk.getUint8(12)],...Array.from(this.raw)])}}st(Br,"type","ihdr"),Pt.set("ihdr",Br),Ct(Nt,"ihdr",[[0,"ImageWidth"],[4,"ImageHeight"],[8,"BitDepth"],[9,"ColorType"],[10,"Compression"],[11,"Filter"],[12,"Interlace"]]),Ct(oe,"ihdr",[[9,{0:"Grayscale",2:"RGB",3:"Palette",4:"Grayscale with Alpha",6:"RGB with Alpha",DEFAULT:"Unknown"}],[10,{0:"Deflate/Inflate",DEFAULT:"Unknown"}],[11,{0:"Adaptive",DEFAULT:"Unknown"}],[12,{0:"Noninterlaced",1:"Adam7 Interlace",DEFAULT:"Unknown"}]]);class qn extends ae{static canHandle(t,e){return t.getUint8(e+1)===226&&t.getUint32(e+4)===1229144927}static findPosition(t,e){let i=super.findPosition(t,e);return i.chunkNumber=t.getUint8(e+16),i.chunkCount=t.getUint8(e+17),i.multiSegment=i.chunkCount>1,i}static handleMultiSegments(t){return function(e){let i=function(r){let s=r[0].constructor,o=0;for(let u of r)o+=u.length;let a=new s(o),l=0;for(let u of r)a.set(u,l),l+=u.length;return a}(e.map(r=>r.chunk.toUint8()));return new Wt(i)}(t)}parse(){return this.raw=new Map,this.parseHeader(),this.parseTags(),this.translate(),this.output}parseHeader(){let{raw:t}=this;this.chunk.byteLength<84&&kt("ICC header is too short");for(let[e,i]of Object.entries(Ga)){e=parseInt(e,10);let r=i(this.chunk,e);r!=="\0\0\0\0"&&t.set(e,r)}}parseTags(){let t,e,i,r,s,{raw:o}=this,a=this.chunk.getUint32(128),l=132,u=this.chunk.byteLength;for(;a--;){if(t=this.chunk.getString(l,4),e=this.chunk.getUint32(l+4),i=this.chunk.getUint32(l+8),r=this.chunk.getString(e,4),e+i>u)return void console.warn("reached the end of the first ICC chunk. Enable options.tiff.multiSegment to read all ICC segments.");s=this.parseTag(r,e,i),s!==void 0&&s!=="\0\0\0\0"&&o.set(t,s),l+=12}}parseTag(t,e,i){switch(t){case"desc":return this.parseDesc(e);case"mluc":return this.parseMluc(e);case"text":return this.parseText(e,i);case"sig ":return this.parseSig(e)}if(!(e+i>this.chunk.byteLength))return this.chunk.getUint8Array(e,i)}parseDesc(t){let e=this.chunk.getUint32(t+8)-1;return Ke(this.chunk.getString(t+12,e))}parseText(t,e){return Ke(this.chunk.getString(t+8,e-8))}parseSig(t){return Ke(this.chunk.getString(t+8,4))}parseMluc(t){let{chunk:e}=this,i=e.getUint32(t+8),r=e.getUint32(t+12),s=t+16,o=[];for(let a=0;a<i;a++){let l=e.getString(s+0,2),u=e.getString(s+2,2),h=e.getUint32(s+4),c=e.getUint32(s+8)+t,m=Ke(e.getUnicodeString(c,h));o.push({lang:l,country:u,text:m}),s+=r}return i===1?o[0].text:o}translateValue(t,e){return typeof t=="string"?e[t]||e[t.toLowerCase()]||t:e[t]||t}}st(qn,"type","icc"),st(qn,"multiSegment",!0),st(qn,"headerLength",18);const Ga={4:we,8:function(n,t){return[n.getUint8(t),n.getUint8(t+1)>>4,n.getUint8(t+1)%16].map(e=>e.toString(10)).join(".")},12:we,16:we,20:we,24:function(n,t){const e=n.getUint16(t),i=n.getUint16(t+2)-1,r=n.getUint16(t+4),s=n.getUint16(t+6),o=n.getUint16(t+8),a=n.getUint16(t+10);return new Date(Date.UTC(e,i,r,s,o,a))},36:we,40:we,48:we,52:we,64:(n,t)=>n.getUint32(t),80:we};function we(n,t){return Ke(n.getString(t,4))}Pt.set("icc",qn),Ct(Nt,"icc",[[4,"ProfileCMMType"],[8,"ProfileVersion"],[12,"ProfileClass"],[16,"ColorSpaceData"],[20,"ProfileConnectionSpace"],[24,"ProfileDateTime"],[36,"ProfileFileSignature"],[40,"PrimaryPlatform"],[44,"CMMFlags"],[48,"DeviceManufacturer"],[52,"DeviceModel"],[56,"DeviceAttributes"],[64,"RenderingIntent"],[68,"ConnectionSpaceIlluminant"],[80,"ProfileCreator"],[84,"ProfileID"],["Header","ProfileHeader"],["MS00","WCSProfiles"],["bTRC","BlueTRC"],["bXYZ","BlueMatrixColumn"],["bfd","UCRBG"],["bkpt","MediaBlackPoint"],["calt","CalibrationDateTime"],["chad","ChromaticAdaptation"],["chrm","Chromaticity"],["ciis","ColorimetricIntentImageState"],["clot","ColorantTableOut"],["clro","ColorantOrder"],["clrt","ColorantTable"],["cprt","ProfileCopyright"],["crdi","CRDInfo"],["desc","ProfileDescription"],["devs","DeviceSettings"],["dmdd","DeviceModelDesc"],["dmnd","DeviceMfgDesc"],["dscm","ProfileDescriptionML"],["fpce","FocalPlaneColorimetryEstimates"],["gTRC","GreenTRC"],["gXYZ","GreenMatrixColumn"],["gamt","Gamut"],["kTRC","GrayTRC"],["lumi","Luminance"],["meas","Measurement"],["meta","Metadata"],["mmod","MakeAndModel"],["ncl2","NamedColor2"],["ncol","NamedColor"],["ndin","NativeDisplayInfo"],["pre0","Preview0"],["pre1","Preview1"],["pre2","Preview2"],["ps2i","PS2RenderingIntent"],["ps2s","PostScript2CSA"],["psd0","PostScript2CRD0"],["psd1","PostScript2CRD1"],["psd2","PostScript2CRD2"],["psd3","PostScript2CRD3"],["pseq","ProfileSequenceDesc"],["psid","ProfileSequenceIdentifier"],["psvm","PS2CRDVMSize"],["rTRC","RedTRC"],["rXYZ","RedMatrixColumn"],["resp","OutputResponse"],["rhoc","ReflectionHardcopyOrigColorimetry"],["rig0","PerceptualRenderingIntentGamut"],["rig2","SaturationRenderingIntentGamut"],["rpoc","ReflectionPrintOutputColorimetry"],["sape","SceneAppearanceEstimates"],["scoe","SceneColorimetryEstimates"],["scrd","ScreeningDesc"],["scrn","Screening"],["targ","CharTarget"],["tech","Technology"],["vcgt","VideoCardGamma"],["view","ViewingConditions"],["vued","ViewingCondDesc"],["wtpt","MediaWhitePoint"]]);const Xn={"4d2p":"Erdt Systems",AAMA:"Aamazing Technologies",ACER:"Acer",ACLT:"Acolyte Color Research",ACTI:"Actix Sytems",ADAR:"Adara Technology",ADBE:"Adobe",ADI:"ADI Systems",AGFA:"Agfa Graphics",ALMD:"Alps Electric",ALPS:"Alps Electric",ALWN:"Alwan Color Expertise",AMTI:"Amiable Technologies",AOC:"AOC International",APAG:"Apago",APPL:"Apple Computer",AST:"AST","AT&T":"AT&T",BAEL:"BARBIERI electronic",BRCO:"Barco NV",BRKP:"Breakpoint",BROT:"Brother",BULL:"Bull",BUS:"Bus Computer Systems","C-IT":"C-Itoh",CAMR:"Intel",CANO:"Canon",CARR:"Carroll Touch",CASI:"Casio",CBUS:"Colorbus PL",CEL:"Crossfield",CELx:"Crossfield",CGS:"CGS Publishing Technologies International",CHM:"Rochester Robotics",CIGL:"Colour Imaging Group, London",CITI:"Citizen",CL00:"Candela",CLIQ:"Color IQ",CMCO:"Chromaco",CMiX:"CHROMiX",COLO:"Colorgraphic Communications",COMP:"Compaq",COMp:"Compeq/Focus Technology",CONR:"Conrac Display Products",CORD:"Cordata Technologies",CPQ:"Compaq",CPRO:"ColorPro",CRN:"Cornerstone",CTX:"CTX International",CVIS:"ColorVision",CWC:"Fujitsu Laboratories",DARI:"Darius Technology",DATA:"Dataproducts",DCP:"Dry Creek Photo",DCRC:"Digital Contents Resource Center, Chung-Ang University",DELL:"Dell Computer",DIC:"Dainippon Ink and Chemicals",DICO:"Diconix",DIGI:"Digital","DL&C":"Digital Light & Color",DPLG:"Doppelganger",DS:"Dainippon Screen",DSOL:"DOOSOL",DUPN:"DuPont",EPSO:"Epson",ESKO:"Esko-Graphics",ETRI:"Electronics and Telecommunications Research Institute",EVER:"Everex Systems",EXAC:"ExactCODE",Eizo:"Eizo",FALC:"Falco Data Products",FF:"Fuji Photo Film",FFEI:"FujiFilm Electronic Imaging",FNRD:"Fnord Software",FORA:"Fora",FORE:"Forefront Technology",FP:"Fujitsu",FPA:"WayTech Development",FUJI:"Fujitsu",FX:"Fuji Xerox",GCC:"GCC Technologies",GGSL:"Global Graphics Software",GMB:"Gretagmacbeth",GMG:"GMG",GOLD:"GoldStar Technology",GOOG:"Google",GPRT:"Giantprint",GTMB:"Gretagmacbeth",GVC:"WayTech Development",GW2K:"Sony",HCI:"HCI",HDM:"Heidelberger Druckmaschinen",HERM:"Hermes",HITA:"Hitachi America",HP:"Hewlett-Packard",HTC:"Hitachi",HiTi:"HiTi Digital",IBM:"IBM",IDNT:"Scitex",IEC:"Hewlett-Packard",IIYA:"Iiyama North America",IKEG:"Ikegami Electronics",IMAG:"Image Systems",IMI:"Ingram Micro",INTC:"Intel",INTL:"N/A (INTL)",INTR:"Intra Electronics",IOCO:"Iocomm International Technology",IPS:"InfoPrint Solutions Company",IRIS:"Scitex",ISL:"Ichikawa Soft Laboratory",ITNL:"N/A (ITNL)",IVM:"IVM",IWAT:"Iwatsu Electric",Idnt:"Scitex",Inca:"Inca Digital Printers",Iris:"Scitex",JPEG:"Joint Photographic Experts Group",JSFT:"Jetsoft Development",JVC:"JVC Information Products",KART:"Scitex",KFC:"KFC Computek Components",KLH:"KLH Computers",KMHD:"Konica Minolta",KNCA:"Konica",KODA:"Kodak",KYOC:"Kyocera",Kart:"Scitex",LCAG:"Leica",LCCD:"Leeds Colour",LDAK:"Left Dakota",LEAD:"Leading Technology",LEXM:"Lexmark International",LINK:"Link Computer",LINO:"Linotronic",LITE:"Lite-On",Leaf:"Leaf",Lino:"Linotronic",MAGC:"Mag Computronic",MAGI:"MAG Innovision",MANN:"Mannesmann",MICN:"Micron Technology",MICR:"Microtek",MICV:"Microvitec",MINO:"Minolta",MITS:"Mitsubishi Electronics America",MITs:"Mitsuba",MNLT:"Minolta",MODG:"Modgraph",MONI:"Monitronix",MONS:"Monaco Systems",MORS:"Morse Technology",MOTI:"Motive Systems",MSFT:"Microsoft",MUTO:"MUTOH INDUSTRIES",Mits:"Mitsubishi Electric",NANA:"NANAO",NEC:"NEC",NEXP:"NexPress Solutions",NISS:"Nissei Sangyo America",NKON:"Nikon",NONE:"none",OCE:"Oce Technologies",OCEC:"OceColor",OKI:"Oki",OKID:"Okidata",OKIP:"Okidata",OLIV:"Olivetti",OLYM:"Olympus",ONYX:"Onyx Graphics",OPTI:"Optiquest",PACK:"Packard Bell",PANA:"Matsushita Electric Industrial",PANT:"Pantone",PBN:"Packard Bell",PFU:"PFU",PHIL:"Philips Consumer Electronics",PNTX:"HOYA",POne:"Phase One A/S",PREM:"Premier Computer Innovations",PRIN:"Princeton Graphic Systems",PRIP:"Princeton Publishing Labs",QLUX:"Hong Kong",QMS:"QMS",QPCD:"QPcard AB",QUAD:"QuadLaser",QUME:"Qume",RADI:"Radius",RDDx:"Integrated Color Solutions",RDG:"Roland DG",REDM:"REDMS Group",RELI:"Relisys",RGMS:"Rolf Gierling Multitools",RICO:"Ricoh",RNLD:"Edmund Ronald",ROYA:"Royal",RPC:"Ricoh Printing Systems",RTL:"Royal Information Electronics",SAMP:"Sampo",SAMS:"Samsung",SANT:"Jaime Santana Pomares",SCIT:"Scitex",SCRN:"Dainippon Screen",SDP:"Scitex",SEC:"Samsung",SEIK:"Seiko Instruments",SEIk:"Seikosha",SGUY:"ScanGuy.com",SHAR:"Sharp Laboratories",SICC:"International Color Consortium",SONY:"Sony",SPCL:"SpectraCal",STAR:"Star",STC:"Sampo Technology",Scit:"Scitex",Sdp:"Scitex",Sony:"Sony",TALO:"Talon Technology",TAND:"Tandy",TATU:"Tatung",TAXA:"TAXAN America",TDS:"Tokyo Denshi Sekei",TECO:"TECO Information Systems",TEGR:"Tegra",TEKT:"Tektronix",TI:"Texas Instruments",TMKR:"TypeMaker",TOSB:"Toshiba",TOSH:"Toshiba",TOTK:"TOTOKU ELECTRIC",TRIU:"Triumph",TSBT:"Toshiba",TTX:"TTX Computer Products",TVM:"TVM Professional Monitor",TW:"TW Casper",ULSX:"Ulead Systems",UNIS:"Unisys",UTZF:"Utz Fehlau & Sohn",VARI:"Varityper",VIEW:"Viewsonic",VISL:"Visual communication",VIVO:"Vivo Mobile Communication",WANG:"Wang",WLBR:"Wilbur Imaging",WTG2:"Ware To Go",WYSE:"WYSE Technology",XERX:"Xerox",XRIT:"X-Rite",ZRAN:"Zoran",Zebr:"Zebra Technologies",appl:"Apple Computer",bICC:"basICColor",berg:"bergdesign",ceyd:"Integrated Color Solutions",clsp:"MacDermid ColorSpan",ds:"Dainippon Screen",dupn:"DuPont",ffei:"FujiFilm Electronic Imaging",flux:"FluxData",iris:"Scitex",kart:"Scitex",lcms:"Little CMS",lino:"Linotronic",none:"none",ob4d:"Erdt Systems",obic:"Medigraph",quby:"Qubyx Sarl",scit:"Scitex",scrn:"Dainippon Screen",sdp:"Scitex",siwi:"SIWI GRAFIKA",yxym:"YxyMaster"},Dr={scnr:"Scanner",mntr:"Monitor",prtr:"Printer",link:"Device Link",abst:"Abstract",spac:"Color Space Conversion Profile",nmcl:"Named Color",cenc:"ColorEncodingSpace profile",mid:"MultiplexIdentification profile",mlnk:"MultiplexLink profile",mvis:"MultiplexVisualization profile",nkpf:"Nikon Input Device Profile (NON-STANDARD!)"};Ct(oe,"icc",[[4,Xn],[12,Dr],[40,Object.assign({},Xn,Dr)],[48,Xn],[80,Xn],[64,{0:"Perceptual",1:"Relative Colorimetric",2:"Saturation",3:"Absolute Colorimetric"}],["tech",{amd:"Active Matrix Display",crt:"Cathode Ray Tube Display",kpcd:"Photo CD",pmd:"Passive Matrix Display",dcam:"Digital Camera",dcpj:"Digital Cinema Projector",dmpc:"Digital Motion Picture Camera",dsub:"Dye Sublimation Printer",epho:"Electrophotographic Printer",esta:"Electrostatic Printer",flex:"Flexography",fprn:"Film Writer",fscn:"Film Scanner",grav:"Gravure",ijet:"Ink Jet Printer",imgs:"Photo Image Setter",mpfr:"Motion Picture Film Recorder",mpfs:"Motion Picture Film Scanner",offs:"Offset Lithography",pjtv:"Projection Television",rpho:"Photographic Paper Printer",rscn:"Reflective Scanner",silk:"Silkscreen",twax:"Thermal Wax Printer",vidc:"Video Camera",vidm:"Video Monitor"}]]);class zn extends ae{static canHandle(t,e,i){return t.getUint8(e+1)===237&&t.getString(e+4,9)==="Photoshop"&&this.containsIptc8bim(t,e,i)!==void 0}static headerLength(t,e,i){let r,s=this.containsIptc8bim(t,e,i);if(s!==void 0)return r=t.getUint8(e+s+7),r%2!=0&&(r+=1),r===0&&(r=4),s+8+r}static containsIptc8bim(t,e,i){for(let r=0;r<i;r++)if(this.isIptcSegmentHead(t,e+r))return r}static isIptcSegmentHead(t,e){return t.getUint8(e)===56&&t.getUint32(e)===943868237&&t.getUint16(e+4)===1028}parse(){let{raw:t}=this,e=this.chunk.byteLength-1,i=!1;for(let r=0;r<e;r++)if(this.chunk.getUint8(r)===28&&this.chunk.getUint8(r+1)===2){i=!0;let s=this.chunk.getUint16(r+3),o=this.chunk.getUint8(r+2),a=this.chunk.getLatin1String(r+5,s);t.set(o,this.pluralizeValue(t.get(o),a)),r+=4+s}else if(i)break;return this.translate(),this.output}pluralizeValue(t,e){return t!==void 0?t instanceof Array?(t.push(e),t):[t,e]:e}}st(zn,"type","iptc"),st(zn,"translateValues",!1),st(zn,"reviveValues",!1),Pt.set("iptc",zn),Ct(Nt,"iptc",[[0,"ApplicationRecordVersion"],[3,"ObjectTypeReference"],[4,"ObjectAttributeReference"],[5,"ObjectName"],[7,"EditStatus"],[8,"EditorialUpdate"],[10,"Urgency"],[12,"SubjectReference"],[15,"Category"],[20,"SupplementalCategories"],[22,"FixtureIdentifier"],[25,"Keywords"],[26,"ContentLocationCode"],[27,"ContentLocationName"],[30,"ReleaseDate"],[35,"ReleaseTime"],[37,"ExpirationDate"],[38,"ExpirationTime"],[40,"SpecialInstructions"],[42,"ActionAdvised"],[45,"ReferenceService"],[47,"ReferenceDate"],[50,"ReferenceNumber"],[55,"DateCreated"],[60,"TimeCreated"],[62,"DigitalCreationDate"],[63,"DigitalCreationTime"],[65,"OriginatingProgram"],[70,"ProgramVersion"],[75,"ObjectCycle"],[80,"Byline"],[85,"BylineTitle"],[90,"City"],[92,"Sublocation"],[95,"State"],[100,"CountryCode"],[101,"Country"],[103,"OriginalTransmissionReference"],[105,"Headline"],[110,"Credit"],[115,"Source"],[116,"CopyrightNotice"],[118,"Contact"],[120,"Caption"],[121,"LocalCaption"],[122,"Writer"],[125,"RasterizedCaption"],[130,"ImageType"],[131,"ImageOrientation"],[135,"LanguageIdentifier"],[150,"AudioType"],[151,"AudioSamplingRate"],[152,"AudioSamplingResolution"],[153,"AudioDuration"],[154,"AudioOutcue"],[184,"JobID"],[185,"MasterDocumentID"],[186,"ShortDocumentID"],[187,"UniqueDocumentID"],[188,"OwnerID"],[200,"ObjectPreviewFileFormat"],[201,"ObjectPreviewFileVersion"],[202,"ObjectPreviewData"],[221,"Prefs"],[225,"ClassifyState"],[228,"SimilarityIndex"],[230,"DocumentNotes"],[231,"DocumentHistory"],[232,"ExifCameraInfo"],[255,"CatalogSets"]]),Ct(oe,"iptc",[[10,{0:"0 (reserved)",1:"1 (most urgent)",2:"2",3:"3",4:"4",5:"5 (normal urgency)",6:"6",7:"7",8:"8 (least urgent)",9:"9 (user-defined priority)"}],[75,{a:"Morning",b:"Both Morning and Evening",p:"Evening"}],[131,{L:"Landscape",P:"Portrait",S:"Square"}]]);let Ui=null;async function Va(){if(!Ui)try{const i=(await import("./joraw-1Lq5hXK7.js")).default;if(typeof i!="function")throw new Error("JoRaw WASM import failed");const r=new URL("/assets/joraw-DraTMNgX.wasm",import.meta.url).href;Ui=i({locateFile:(s,o)=>s.endsWith("joraw.wasm")?r:o+s})}catch(e){throw console.error("Failed to load joraw.js:",e),e}const n=await Ui,t=n.LibRaw||n.JoRaw;if(!t)throw new Error("JoRaw class not found");return t}const Xa=async n=>{var i,r,s,o,a,l,u,h,c;const t=await Va(),e=new t;try{if(await e.open(n,{}),typeof e.getRawImage!="function")throw new Error("WASM mismatch");const m=e.getRawImage();let d=new Uint16Array(m.data);const f=await e.metadata(!0);let p={...f};try{const S=await Ba.parse(n.buffer);S&&(p={...p,...S})}catch(S){console.warn("exifr parsing failed for RAW buffer",S)}const g=((i=f.idata)==null?void 0:i.filters)||0,y=((r=f.idata)==null?void 0:r.colors)||0,x=g===0&&y===3,b=g===9;let _=[0,0,0,0],P=!1;if(e.getBlackLevels)try{const S=e.getBlackLevels();S.dng_cblack&&S.dng_cblack.length===4&&Array.from(S.dng_cblack).some(w=>w>0)?(_=Array.from(S.dng_cblack).map(Number),P=!0):S.cblack&&S.cblack.length===4&&Array.from(S.cblack).some(w=>w>0)?(_=Array.from(S.cblack).map(Number),P=!0):typeof S.black=="number"&&S.black>0&&(_=[S.black,S.black,S.black,S.black],P=!0)}catch(S){console.warn("getBlackLevels binding failed",S)}if(!P){let S=[];if((s=f.color_data)!=null&&s.cblack_rawpy_style)S=f.color_data.cblack_rawpy_style;else if((a=(o=f.color_data)==null?void 0:o.dng_levels)!=null&&a.dng_cblack)S=f.color_data.dng_levels.dng_cblack;else if(((l=f.black_level_per_channel)==null?void 0:l.length)>=4)S=f.black_level_per_channel;else if(((u=f.cblack)==null?void 0:u.length)>=4)S=f.cblack;else if(((c=(h=f.color)==null?void 0:h.cblack)==null?void 0:c.length)>=4)S=f.color.cblack;else{const w=f.black_level||f.color_data&&f.color_data.black||0;S=[w,w,w,w]}_=[Number(S[0])||0,Number(S[1])||0,Number(S[2])||0,Number(S[3])||0]}return{data:d,width:m.width,height:m.height,bayerPattern:f.color_desc||"RGGB",blackLevels:_,whiteLevel:f.white_level||16383,metadata:p,isThreePlane:x,threePlaneTransfer:x?"linear":void 0,isXTrans:b}}finally{e.delete?e.delete():e.close()}};async function za(n){if(ea(n)){const i=await oa(n);if(!i)throw new Error("Sony cRAW HQ decoder did not return image data.");return i.rawImageData}if(pa(n)){const i=await ma(n);if(!i)throw new Error("Nikon HE decoder did not return image data.");return i}const e=new Uint8Array(n);return Xa(e)}function Ya(n,t,e){const i=fs(n,e),r=Math.floor(t.x),s=Math.floor(t.y),o=Math.floor(t.w),a=Math.floor(t.h),l=new Uint16Array(o*a);for(let u=0;u<a;u++){const h=s+u,c=u*o;for(let m=0;m<o;m++)l[c+m]=i(r+m,h)}return{data:l,width:o,height:a}}function Wa(n,t,e){const i=Ha(n,t,e);if(!i)return null;const r=n.width,s=n.height,o=new Uint16Array(r*s);for(let a=0;a<s;a++){const l=a*r;for(let u=0;u<r;u++)o[l+u]=i(u,a)}return{kind:"u16-mono",data:o,width:r,height:s}}function Ha(n,t,e){return t.renderMode==="advanced-zero-dep"&&t.advancedZeroDep?fs(n,t,e):t.renderMode==="zero-dependency"?ja(n,t,e):null}function fs(n,t,e){if(!t.advancedZeroDep)throw new Error("Unmixing settings not found in DisplaySettings.");const{bg:i,fg:r}=t.advancedZeroDep,s=ps(e,t.advancedZeroDep.bl),{data:o,width:a,whiteLevel:l}=n,u=i.map((f,p)=>Math.max(0,f-s[p])),h=r.map((f,p)=>Math.max(0,f-s[p])),c=(u[1]+u[3])/2,m=(h[1]+h[3])/2,d=Math.pow(2,t.exposure);return(f,p)=>{if(f<0||p<0||f>=a||p>=n.height)return 0;const g=p%2,y=f%2;let x=0;!g&&!y?x=0:!g&&y?x=1:g&&!y?x=3:x=2;const b=o[p*a+f],_=s[x],P=u[x],S=h[x],w=Math.max(b-_,0),C=S-P||1e-9,A=(w-P)/C;let k;return A<0?k=w*(c/Math.max(P,1e-9)):A>1?k=w*(m/Math.max(S,1e-9)):k=(1-A)*c+A*m,k*=d,Math.max(0,Math.min(65535,Math.round(k)))}}function ja(n,t,e){const{data:i,width:r,height:s,whiteLevel:o}=n,a=ps(e,t.blackLevel||[0,0,0,0]),l=qa(n.bayerPattern),u=t.wbGains?t.wbGains[0]:1,h=t.wbGains?t.wbGains[1]:1,c=Math.pow(2,t.exposure||0);return(m,d)=>{if(m<0||d<0||m>=r||d>=s)return 0;const f=ms(m,d),p=Ka(l,m,d),g=i[d*r+m],y=a[f];let x=(g-y)/Math.max(1,o-y);return x=Math.max(0,Math.min(1,x)),x*=c,p==="R"?x*=u:p==="B"&&(x*=h),Qa(x)}}function ps(n,t){if(typeof n=="number"&&Number.isFinite(n)){const e=Math.max(0,n);return[e,e,e,e]}return Array.isArray(n)&&n.length===4?[Math.max(0,n[0]??0),Math.max(0,n[1]??0),Math.max(0,n[2]??0),Math.max(0,n[3]??0)]:[Math.max(0,t[0]??0),Math.max(0,t[1]??0),Math.max(0,t[2]??0),Math.max(0,t[3]??0)]}function Qa(n){return Math.max(0,Math.min(65535,Math.round(Math.max(0,Math.min(1,n))*65535)))}function qa(n){const t=(n||"RGGB").toUpperCase().trim();return t.length>=4&&/^[RGB]{4}$/.test(t.slice(0,4))?t.slice(0,4):"RGGB"}function ms(n,t){return(t&1)<<1|n&1}function Ka(n,t,e){return n[ms(t,e)]}const ft=(n,t=0)=>({real:n,imag:t}),fn=(n,t)=>({real:n.real+t.real,imag:n.imag+t.imag}),hn=(n,t)=>({real:n.real-t.real,imag:n.imag-t.imag}),It=(n,t)=>({real:n.real*t.real-n.imag*t.imag,imag:n.real*t.imag+n.imag*t.real}),Le=(n,t)=>{const e=t.real*t.real+t.imag*t.imag;return e===0?ft(0):{real:(n.real*t.real+n.imag*t.imag)/e,imag:(n.imag*t.real-n.real*t.imag)/e}},qe=n=>Math.hypot(n.real,n.imag),gs=n=>{const t=qe(n);if(t===0)return ft(0);const e=Math.sqrt(t),i=Math.atan2(n.imag,n.real);return ft(e*Math.cos(i/2),e*Math.sin(i/2))};function $a(n,t){const e=n.length-1;if(e<0)return{p:ft(0),dp:ft(0),d2p:ft(0)};let i=ft(n[e].real,n[e].imag),r=ft(0),s=ft(0);for(let o=e-1;o>=0;o--)s=fn(It(r,ft(2)),It(t,s)),r=fn(i,It(t,r)),i=fn(ft(n[o].real,n[o].imag),It(t,i));return{p:i,dp:r,d2p:s}}function Wi(n,t,e=80){const r=n.length-1;if(r<=0)return{root:t,iterations:0};if(r===1)return{root:Le(It(n[0],ft(-1)),n[1]),iterations:0};let s=ft(t.real,t.imag);for(let o=0;o<e;o++){const{p:a,dp:l,d2p:u}=$a(n,s);if(qe(a)<1e-14)return{root:s,iterations:o};const h=Le(l,a),c=It(h,h),m=hn(c,Le(u,a)),d=ft(r),f=ft(r-1),p=hn(It(d,m),It(h,h)),g=gs(It(f,p)),y=fn(h,g),x=hn(h,g),b=qe(y)>qe(x)?y:x;if(qe(b)<1e-14)return{root:s,iterations:o};const _=Le(d,b),P=hn(s,_);if(qe(_)<1e-14*qe(P))return{root:P,iterations:o+1};s=P}return{root:s,iterations:e}}function Ja(n,t){const e=n.length-1;if(e<=0)return[ft(0)];if(e===1)return[ft(n[0].real,n[0].imag)];const i=new Array(e);i[e-1]=ft(n[e].real,n[e].imag);for(let r=e-2;r>=0;r--){const s=ft(n[r+1].real,n[r+1].imag),o=i[r+1];i[r]=fn(s,It(t,o))}return i}function Za(n){const t=n.length-1;if(t<=0)return[];if(t===1)return[Le(It(n[0],ft(-1)),n[1])];const e=[];let i=n.map(s=>ft(s.real,s.imag)),r=t*5;for(;i.length>2&&r-- >0;){const s=ft(.3+Math.random()*.7,.3+Math.random()*.7),{root:o}=Wi(i,s,100),a=Wi(n,o,20);e.push(a.root);const l=Ja(i,o);if(l.length>=i.length){console.warn("polyDeflate did not reduce degree, breaking");break}i=l}if(i.length===2)e.push(Le(It(i[0],ft(-1)),i[1]));else if(i.length===3){const s=i[2],o=i[1],a=i[0],l=hn(It(o,o),It(It(ft(4),s),a)),u=gs(l),h=It(ft(2),s),c=Le(hn(It(ft(-1),o),u),h),m=Le(fn(It(ft(-1),o),u),h);e.push(c,m)}return e}function to(n,t,e){const i=[ft(n),ft(-1),ft(n*t),ft(0),ft(n*e)],r=Za(i);if(r.length===0)return console.warn("laguerreSmallestPositiveRoot: no roots found"),n;let s=1/0,o=!1;for(const l of r)Math.abs(l.imag)<1e-10&&l.real>0&&l.real<s&&(s=l.real,o=!0);return o?Wi(i,ft(s,0),20).root.real:(console.warn("laguerreSmallestPositiveRoot: no positive real root found"),n)}function eo(n,t,e){if(Math.abs(t)<1e-10&&Math.abs(e)<1e-10)return n;if(n<1e-10)return 0;if(Math.abs(e)<1e-10){const i=-1/(t*n),r=1/t,s=i*i-4*r;if(s<0)return n;const o=Math.sqrt(s),a=-.5*(i+Math.sign(i)*o),l=a,u=r/a;return l>0&&u>0?Math.min(l,u):l>0?l:u>0?u:n}try{return to(n,t,e)}catch(i){return console.error("laguerreSmallestPositiveRoot failed:",i),n}}function no(n,t,e,i,r){const s=n.x-t.x,o=n.y-t.y,a=Math.hypot(s,o)/Math.max(1e-12,e);if(a<1e-12)return{x:n.x,y:n.y};const l=a*a,u=1+(i+r*l)*l;return{x:s/u+t.x,y:o/u+t.y}}function io(n,t,e,i,r){const s=n.x-t.x,o=n.y-t.y,a=Math.hypot(s,o)/Math.max(1e-12,e);if(a<1e-12)return{x:t.x,y:t.y};const u=eo(a,i,r)/a;return{x:t.x+s*u,y:t.y+o*u}}function Ut(n,t){const e=no(n,{x:t.principalX,y:t.principalY},t.radiusNorm,t.k1,t.k2);return{x:e.x+(t.correctedOffsetX??0),y:e.y+(t.correctedOffsetY??0)}}function Pe(n,t){const e={x:n.x-(t.correctedOffsetX??0),y:n.y-(t.correctedOffsetY??0)};return io(e,{x:t.principalX,y:t.principalY},t.radiusNorm,t.k1,t.k2)}const ro=`
attribute vec2 a_position;

void main() {
  gl_Position = vec4(a_position, 0.0, 1.0);
}
`,so=`
precision highp float;
precision mediump int;

uniform sampler2D u_source;
uniform vec2 u_size;

vec2 pixelUv(float x, float y) {
  return (vec2(x, y) + 0.5) / u_size;
}

float sampleGray(float x, float y) {
  return floor(texture2D(u_source, pixelUv(x, y)).r * 255.0 + 0.5);
}

void main() {
  float x = floor(gl_FragCoord.x - 0.5);
  float y = floor(gl_FragCoord.y - 0.5);
  float sum = 0.0;
  float count = 0.0;

  for (int oy = -1; oy <= 1; oy++) {
    for (int ox = -1; ox <= 1; ox++) {
      float sx = x + float(ox);
      float sy = y + float(oy);
      if (sx < 0.0 || sy < 0.0 || sx >= u_size.x || sy >= u_size.y) continue;
      sum += sampleGray(sx, sy);
      count += 1.0;
    }
  }

  float blurred = floor(sum / max(count, 1.0) + 0.5);
  float encoded = blurred / 255.0;
  gl_FragColor = vec4(encoded, encoded, encoded, 1.0);
}
`,ao=`
precision highp float;
precision mediump int;

uniform sampler2D u_blurred;
uniform vec2 u_size;

vec2 pixelUv(float x, float y) {
  return (vec2(x, y) + 0.5) / u_size;
}

float sampleBlurred(float x, float y) {
  return floor(texture2D(u_blurred, pixelUv(x, y)).r * 255.0 + 0.5);
}

vec2 packSigned16(float value) {
  float shifted = floor(clamp(value + 32768.0, 0.0, 65535.0) + 0.5);
  float lo = mod(shifted, 256.0);
  float hi = floor(shifted / 256.0);
  return vec2(lo, hi) / 255.0;
}

void main() {
  float x = floor(gl_FragCoord.x - 0.5);
  float y = floor(gl_FragCoord.y - 0.5);

  if (x < 1.0 || y < 1.0 || x >= u_size.x - 1.0 || y >= u_size.y - 1.0) {
    vec2 zero = packSigned16(0.0);
    gl_FragColor = vec4(zero, zero);
    return;
  }

  float tl = sampleBlurred(x - 1.0, y - 1.0);
  float tc = sampleBlurred(x, y - 1.0);
  float tr = sampleBlurred(x + 1.0, y - 1.0);
  float ml = sampleBlurred(x - 1.0, y);
  float mr = sampleBlurred(x + 1.0, y);
  float bl = sampleBlurred(x - 1.0, y + 1.0);
  float bc = sampleBlurred(x, y + 1.0);
  float br = sampleBlurred(x + 1.0, y + 1.0);

  float dx = (-tl - 2.0 * ml - bl) + (tr + 2.0 * mr + br);
  float dy = (-tl - 2.0 * tc - tr) + (bl + 2.0 * bc + br);

  vec2 packedDx = packSigned16(dx);
  vec2 packedDy = packSigned16(dy);
  gl_FragColor = vec4(packedDx, packedDy);
}
`;class oo{constructor(){bt(this,"canvas",null);bt(this,"gl",null);bt(this,"blurProgram",null);bt(this,"sobelProgram",null);bt(this,"positionBuffer",null);bt(this,"blurUniforms",null);bt(this,"sobelUniforms",null);bt(this,"resources",null);bt(this,"initialized",!1);bt(this,"unavailable",!1);bt(this,"maxTextureSize",0)}compute(t,e,i){if(!this.initialized&&!this.init())return null;const r=this.gl,s=this.blurProgram,o=this.sobelProgram,a=this.blurUniforms,l=this.sobelUniforms;if(!r||!s||!o||!a||!l||!this.positionBuffer||!this.canvas||e<=2||i<=2||e>this.maxTextureSize||i>this.maxTextureSize)return null;const u=this.ensureResources(e,i);if(!u)return null;this.canvas.width=e,this.canvas.height=i,r.viewport(0,0,e,i),r.disable(r.BLEND),r.pixelStorei(r.UNPACK_ALIGNMENT,1),r.pixelStorei(r.PACK_ALIGNMENT,1),r.bindBuffer(r.ARRAY_BUFFER,this.positionBuffer),r.activeTexture(r.TEXTURE0),r.bindTexture(r.TEXTURE_2D,u.sourceTexture),r.texImage2D(r.TEXTURE_2D,0,r.LUMINANCE,e,i,0,r.LUMINANCE,r.UNSIGNED_BYTE,t),r.useProgram(s),r.enableVertexAttribArray(0),r.vertexAttribPointer(0,2,r.FLOAT,!1,0,0),r.uniform2f(a.size,e,i),r.uniform1i(a.source,0),r.bindFramebuffer(r.FRAMEBUFFER,u.blurFramebuffer),r.bindTexture(r.TEXTURE_2D,u.sourceTexture),r.drawArrays(r.TRIANGLES,0,6);const h=new Uint8Array(e*i*4);r.readPixels(0,0,e,i,r.RGBA,r.UNSIGNED_BYTE,h),r.useProgram(o),r.uniform2f(l.size,e,i),r.uniform1i(l.blurred,0),r.bindFramebuffer(r.FRAMEBUFFER,u.sobelFramebuffer),r.bindTexture(r.TEXTURE_2D,u.blurTexture),r.drawArrays(r.TRIANGLES,0,6);const c=new Uint8Array(e*i*4);r.readPixels(0,0,e,i,r.RGBA,r.UNSIGNED_BYTE,c),r.disableVertexAttribArray(0),r.bindFramebuffer(r.FRAMEBUFFER,null),r.bindBuffer(r.ARRAY_BUFFER,null),r.bindTexture(r.TEXTURE_2D,null);const m=new Uint8Array(e*i);for(let g=0,y=0;g<m.length;g++,y+=4)m[g]=h[y];const d=new Float32Array(e*i),f=new Float32Array(e*i),p=new Float32Array(e*i);for(let g=0,y=0;g<d.length;g++,y+=4){const x=(c[y]|c[y+1]<<8)-32768,b=(c[y+2]|c[y+3]<<8)-32768;d[g]=x,f[g]=b,p[g]=Math.sqrt(x*x+b*b)}return{blurredGray:m,gx:d,gy:f,magnitude:p}}init(){if(this.initialized&&this.gl&&this.blurProgram&&this.sobelProgram)return!0;if(this.unavailable)return!1;const t=this.createCanvas();if(!t)return this.unavailable=!0,!1;const e=t.getContext("webgl",{alpha:!1,antialias:!1,depth:!1,stencil:!1,premultipliedAlpha:!1,preserveDrawingBuffer:!1});if(!e)return this.unavailable=!0,!1;const i=this.compileShader(e,e.VERTEX_SHADER,ro),r=this.compileShader(e,e.FRAGMENT_SHADER,so),s=this.compileShader(e,e.FRAGMENT_SHADER,ao);if(!i||!r||!s)return i&&e.deleteShader(i),r&&e.deleteShader(r),s&&e.deleteShader(s),this.unavailable=!0,!1;const o=this.createProgram(e,i,r),a=this.createProgram(e,i,s);if(e.deleteShader(i),e.deleteShader(r),e.deleteShader(s),!o||!a)return o&&e.deleteProgram(o),a&&e.deleteProgram(a),this.unavailable=!0,!1;const l=e.createBuffer();return l?(e.bindBuffer(e.ARRAY_BUFFER,l),e.bufferData(e.ARRAY_BUFFER,new Float32Array([-1,-1,1,-1,-1,1,-1,1,1,-1,1,1]),e.STATIC_DRAW),e.bindBuffer(e.ARRAY_BUFFER,null),this.canvas=t,this.gl=e,this.blurProgram=o,this.sobelProgram=a,this.positionBuffer=l,this.blurUniforms={source:e.getUniformLocation(o,"u_source"),size:e.getUniformLocation(o,"u_size")},this.sobelUniforms={blurred:e.getUniformLocation(a,"u_blurred"),size:e.getUniformLocation(a,"u_size")},this.maxTextureSize=Number(e.getParameter(e.MAX_TEXTURE_SIZE)||0),this.initialized=!0,!0):(e.deleteProgram(o),e.deleteProgram(a),this.unavailable=!0,!1)}createCanvas(){return typeof OffscreenCanvas<"u"?new OffscreenCanvas(1,1):typeof document<"u"?document.createElement("canvas"):null}ensureResources(t,e){const i=this.gl;if(!i)return null;if(this.resources&&this.resources.width===t&&this.resources.height===e)return this.resources;this.disposeResources();const r=this.createTexture(i.LUMINANCE,t,e,i.LUMINANCE,i.UNSIGNED_BYTE,null),s=this.createTexture(i.RGBA,t,e,i.RGBA,i.UNSIGNED_BYTE,null),o=this.createTexture(i.RGBA,t,e,i.RGBA,i.UNSIGNED_BYTE,null),a=this.createFramebuffer(s),l=this.createFramebuffer(o);return!r||!s||!o||!a||!l?(r&&i.deleteTexture(r),s&&i.deleteTexture(s),o&&i.deleteTexture(o),a&&i.deleteFramebuffer(a),l&&i.deleteFramebuffer(l),null):(this.resources={width:t,height:e,sourceTexture:r,blurTexture:s,sobelTexture:o,blurFramebuffer:a,sobelFramebuffer:l},this.resources)}createTexture(t,e,i,r,s,o){const a=this.gl;if(!a)return null;const l=a.createTexture();return l?(a.bindTexture(a.TEXTURE_2D,l),a.texParameteri(a.TEXTURE_2D,a.TEXTURE_MIN_FILTER,a.NEAREST),a.texParameteri(a.TEXTURE_2D,a.TEXTURE_MAG_FILTER,a.NEAREST),a.texParameteri(a.TEXTURE_2D,a.TEXTURE_WRAP_S,a.CLAMP_TO_EDGE),a.texParameteri(a.TEXTURE_2D,a.TEXTURE_WRAP_T,a.CLAMP_TO_EDGE),a.texImage2D(a.TEXTURE_2D,0,t,e,i,0,r,s,o),a.bindTexture(a.TEXTURE_2D,null),l):null}createFramebuffer(t){const e=this.gl;if(!e||!t)return null;const i=e.createFramebuffer();if(!i)return null;e.bindFramebuffer(e.FRAMEBUFFER,i),e.framebufferTexture2D(e.FRAMEBUFFER,e.COLOR_ATTACHMENT0,e.TEXTURE_2D,t,0);const r=e.checkFramebufferStatus(e.FRAMEBUFFER);return e.bindFramebuffer(e.FRAMEBUFFER,null),r!==e.FRAMEBUFFER_COMPLETE?(e.deleteFramebuffer(i),null):i}compileShader(t,e,i){const r=t.createShader(e);return r?(t.shaderSource(r,i),t.compileShader(r),t.getShaderParameter(r,t.COMPILE_STATUS)?r:(console.error("[SFR Auto Detect WebGL] shader compile failed",t.getShaderInfoLog(r)),t.deleteShader(r),null)):null}createProgram(t,e,i){const r=t.createProgram();return r?(t.attachShader(r,e),t.attachShader(r,i),t.bindAttribLocation(r,0,"a_position"),t.linkProgram(r),t.getProgramParameter(r,t.LINK_STATUS)?r:(console.error("[SFR Auto Detect WebGL] program link failed",t.getProgramInfoLog(r)),t.deleteProgram(r),null)):null}disposeResources(){const t=this.gl,e=this.resources;if(!t||!e){this.resources=null;return}t.deleteTexture(e.sourceTexture),t.deleteTexture(e.blurTexture),t.deleteTexture(e.sobelTexture),t.deleteFramebuffer(e.blurFramebuffer),t.deleteFramebuffer(e.sobelFramebuffer),this.resources=null}}const lo=new oo,Hi={gradientPercentiles:[.82,.88,.92,.95,.98,.995],downsampleMaxSide:1600,minComponentAreaRatio:15e-6,maxComponentAreaRatio:.35,minComponentAreaPx:20,minEdgePoints:24,extentQuantileLow:.02,extentQuantileHigh:.98,cornerTrimRatio:.18,minSpanPx:8,maxAspectRatio:2,bandScale:.16,bandMinPx:1.75,bandMaxPx:14,minPointContrast:6,minSidePoints:3,minCoverageRatio:.15,minCenterCoverageRatio:.2,filterBlockPurity:!0,innerPurityStdScale:1.5,outerMeanSpreadLimit:51,minAxisDot:.6,residualLimitFloor:.01,residualLimitScale:.25,minQuadArea:48,minSideLength:10,minOuterContrast:5,sampleHalfWidthRatio:.25};function co(n,t,e,i,r,s){const o=n.width,a=n.height,l=(n.bayerPattern||"RGGB").toUpperCase(),u=[],h=[],c=s!=null&&s.correctedRect?Lt*2:Lt,m=Math.max(1,Math.min(r,c)),d=e.p2.x-e.p1.x,f=e.p2.y-e.p1.y,p=Math.hypot(d,f);if(!Number.isFinite(p)||p<=1e-6)return null;const g=d/p,y=f/p,x=-y,b=g,_={x:(e.p1.x+e.p2.x)*.5,y:(e.p1.y+e.p2.y)*.5},P=s!=null&&s.correctedRect?Mt(s.correctedRect,o,a):Mt(Gt(xe(e,c*4+2)??[e.p1,e.p2],2),o,a);if(!P)return null;const S=(s==null?void 0:s.correctedScanlinesOverride)??(s!=null&&s.distortedRect?or(Mt(s.distortedRect,o,a)??s.distortedRect,t,o,a):Ns(P,e,Math.max(1,i),m*4+.5,o,a));if(!S||S.size===0)return null;const w=Rs(S,t,o,a);if(w.size===0)return null;const C=!lr(t);for(const[k,F]of w){if(k<0||k>=a)continue;const M=k*o;for(let v=F.start;v<=F.end;v++){if(v<0||v>=o||!vt(v,k,l,s==null?void 0:s.greenPhase))continue;const I={x:v,y:k},R=Ut(I,t);if(!Number.isFinite(R.x)||!Number.isFinite(R.y)||Math.round(R.x)<0||Math.round(R.x)>=o||Math.round(R.y)<0||Math.round(R.y)>=a)continue;const U=R.x-_.x,Y=R.y-_.y,B=U*g+Y*y;let L=U*x+Y*b;if(C){const O=Ls(B,g,y,_,I,t);if(!O)continue;const D=.5*(O.a+O.b),X=qt(O.a,g,y,_,t),et=qt(D,g,y,_,t),V=qt(O.b,g,y,_,t),$=ur({x:O.a,y:Math.hypot(X.x-I.x,X.y-I.y)},{x:D,y:Math.hypot(et.x-I.x,et.y-I.y)},{x:O.b,y:Math.hypot(V.x-I.x,V.y-I.y)});if(!Number.isFinite($))continue;const H=cr($,g,y,_,t),it=Math.hypot(H.x,H.y);if(!Number.isFinite(it)||it<=1e-9)continue;const at=H.x/it,J=-(H.y/it),Q=at,T=qt($,g,y,_,t);L=(I.x-T.x)*J+(I.y-T.y)*Q}!Number.isFinite(B)||Math.abs(B)>Math.max(1,i)||!Number.isFinite(L)||Math.abs(L)>m||(u.push(L),h.push(Math.max(0,n.data[M+v]-xn(s==null?void 0:s.blackLevel,v,k))))}}if(u.length<8)return null;const A=Math.abs(d)>=Math.abs(f)?1:2;return s!=null&&s.forceLegacyModel?An(u,h,A,c):bn(u,h,A,c)}function nr(n,t){const e=n.length;let i=0,r=0,s=0,o=0;for(let l=0;l<e;l++)i+=n[l],r+=t[l],s+=n[l]*t[l],o+=n[l]*n[l];const a=e*o-i*i;return a===0?{slope:0,intercept:0}:{slope:(e*s-i*r)/a,intercept:(r*o-i*s)/a}}function uo(n,t){const e=n.length,i=new Array(e).fill(0),r=2*t;for(let s=0;s<e;s++){const o=s>0?n[s-1]:n[0],a=s<e-1?n[s+1]:n[e-1];i[s]=(a-o)/r}return i}const Zt=-1e7,pn=13,se=512,Et=8,ui=1/Et,Lt=28,ho=[[0,0,0,0,0,-.085714285714286,.342857142857143,.485714285714286,.342857142857143,-.085714285714286,0,0,0,0,0],[0,0,0,0,-.095238095238095,.142857142857143,.285714285714286,.333333333333333,.285714285714286,.142857142857143,-.095238095238095,0,0,0,0],[0,0,0,-.090909090909091,.060606060606061,.168831168831169,.233766233766234,.255411255411255,.233766233766234,.168831168831169,.060606060606061,-.090909090909091,0,0,0],[0,0,-.083916083916084,.020979020979021,.102564102564103,.160839160839161,.195804195804196,.207459207459208,.195804195804196,.160839160839161,.102564102564103,.020979020979021,-.083916083916084,0,0],[0,-.076923076923077,0,.062937062937063,.111888111888112,.146853146853147,.167832167832168,.174825174825175,.167832167832168,.146853146853147,.111888111888112,.062937062937063,0,-.076923076923077,0],[-.070588235294118,-.011764705882353,.038009049773756,.078733031674208,.110407239819004,.133031674208145,.146606334841629,.151131221719457,.146606334841629,.133031674208145,.110407239819004,.078733031674208,.038009049773756,-.011764705882353,-.070588235294118]];function mn(n,t,e,i=1){const r=Math.max(1e-6,e*.5),s=Math.max(1e-6,t*i),o=Math.exp(-s*r),a=1-o;if(!Number.isFinite(a)||Math.abs(a)<=1e-9)return Math.abs(n)<=r?1:0;if(Math.abs(n)<r){const l=2-2*o*Math.cosh(s*n),u=2*Math.sinh(s*r)*a;return!Number.isFinite(u)||Math.abs(u)<=1e-9?0:l/u}return Math.exp(-s*Math.abs(n))/a}function ys(n,t,e,i,r,s){const o=n.length;if(o===0)return[];if(s<1)return n;const a=Math.min(s,32),l=new Array(o).fill(0);l[0]=n[0];for(let d=1;d<o;d++)l[d]=l[d-1]+n[d];const u=(d,f)=>{const p=Math.max(0,d),g=Math.min(o-1,f);return g<p?n[Math.max(0,Math.min(o-1,d))]??0:(l[g]-(p>0?l[p-1]:0))/(g-p+1)},h=a*2,c=a,m=1;for(let d=Math.max(t+c,i-h);d<i;d++){const f=Math.max(m,Math.trunc((i-d)*c/Math.max(1,h)));n[d]=u(d-f,d+f)}for(let d=Math.min(r+h-1,e-c-1);d>r;d--){const f=Math.max(m,Math.trunc((d-r)*c/Math.max(1,h)));n[d]=u(d-f,d+f)}for(let d=c+1;d<i-h;d++)n[d]=u(d-c,d+c);for(let d=Math.min(r+h,e-c-1);d<o-c-1;d++)n[d]=u(d-c,d+c);return n}function xs(n){return!Number.isFinite(n)||Math.abs(n)<=1e-9?1:Math.sin(n)/n}let Yn=null,Bi=null;function fo(){if(Yn)return Yn;const n=.625,t=1/128,e=Math.max(16,Math.round(n*2/t)+1),i=[],r=[];for(let c=0;c<e;c++){const m=-n+c*t;i.push(m),r.push(Math.abs(m)<=n?mn(m,pn,.125,1):0)}const s=4,o=1/1024,a=Math.round(s/o)+1,l=new Array(a).fill(0),u=new Array(a).fill(1);let h=0;for(let c=0;c<r.length;c++)h+=r[c];h=Math.max(1e-9,h);for(let c=0;c<a;c++){const m=c*o;l[c]=m;let d=0;for(let f=0;f<i.length;f++)d+=r[f]*Math.cos(2*Math.PI*m*i[f]);u[c]=Math.max(1e-6,Math.abs(d)/h)}return Yn={freqs:l,values:u},Yn}function po(n,t){const e=Math.max(1e-6,xs(Math.PI*n*t)),i=fo(),r=Ze(Math.max(0,Math.min(i.freqs[i.freqs.length-1],n)),i.freqs,i.values);return Math.max(1e-6,e*r)}function mo(){if(Bi)return Bi;const n=new Array(se/16*4).fill(1),t=se*16,e=new Float32Array(t);for(let s=0;s<t;s++){const o=(s-t/2)/(16*Et);e[s]=Math.abs(o)<=.625?mn(o,pn,ui,1):0}const i=new Fn(t);i.transform(e);const r=Math.max(1e-9,Math.abs(i._real[0]));n[0]=1;for(let s=1;s<n.length;s++){const o=xs(Math.PI*s/256),a=Math.max(1e-6,Math.hypot(i._real[s],i._imag[s])/r);n[s]=Math.max(1e-6,o*a)}return Bi=n,n}function go(n,t,e){if(n.length===0||t.length!==n.length||!(e>0))return null;const i=n[0],r=n[n.length-1],s=Math.floor((r-i)/e);if(s<16)return null;const o=Math.max(0,Math.min(s-1,Math.round(-i/e))),a=Math.max(2,Math.round(2/e)),l=Math.max(1,Math.round(.5/e)),u=5,h=.125/Math.max(e,1e-6),c=E=>i+E*e,d=((E,J)=>{let Q=Math.max(0,E),T=Math.min(s-1,J);if(T-Q<8)return null;let j=0;for(;;){const G=new Array(s).fill(0),z=new Array(s).fill(0);for(let K=0;K<n.length;K++){const Ft=Math.max(0,Math.min(s-1,Math.trunc((n[K]-i)/e))),Tt=Math.max(Q,Ft-u),xt=Math.min(T-1,Ft+u);for(let wt=Tt;wt<=xt;wt++){const ee=c(wt),$t=Math.max(0,1-Math.abs((n[K]-ee)*1.75*h));$t<=0||(G[wt]+=t[K]*$t,z[wt]+=$t)}}const Z=new Array(s).fill(Zt);let tt=0,nt=0,ct=0,dt=0,rt=-1,Bt=-1;const Vt=Math.max(o-Math.round(s/8),Q+a),W=Math.min(o+Math.round(s/8),T-a);for(let K=Math.max(0,Q-1);K<=Math.min(s-1,T+1);K++)z[K]>0&&(Z[K]=G[K]/z[K],K<Vt&&(tt+=Z[K],ct++),K>W&&(nt+=Z[K],dt++),rt<0&&(rt=K),Bt=K);if(rt<0||Bt<0||ct===0||dt===0)return null;for(let K=rt-1;K>=0;K--)Z[K]=Z[rt];for(let K=Bt+1;K<s;K++)Z[K]=Z[Bt];const ut=Math.max(2,a),ot=new Array(s).fill(0);let pt=o;for(let K=ut+1;K<s-1-ut;K++){let Ft=0,Tt=0;for(let xt=-ut;xt<=ut;xt++)Tt+=Z[K+xt]*xt,Ft+=xt*xt;ot[K]=Ft>0?Tt/Ft:0,Math.abs(ot[K])>Math.abs(ot[pt]??0)&&K>Q+ut&&K<T-ut-1&&(pt=K)}const mt=Math.max(1,Math.round(2/e)),te=Math.max(mt+1,Math.round(12/e)),Dt=Math.abs(pt-o);if(Dt>mt&&Dt<te)return null;let _t=0;for(let K=Math.max(0,o-ut);K<=Math.min(s-1,o+ut);K++)Math.abs(ot[K])>Math.abs(_t)&&(_t=ot[K]);if(!Number.isFinite(_t)||Math.abs(_t)<=1e-9)return null;const le=Math.abs(_t*.001);let Kt=!1,ke=!1;for(let K=Math.max(0,o-a);K>=Q+l;K--)if(ot[K]*_t<0&&Math.abs(ot[K])>le){let Ft=0,Tt=0,xt=0;for(let wt=K;wt>=Q;wt--)ot[wt]*_t<0&&(Ft++,Tt=Math.max(Tt,Math.abs(ot[wt]))),xt++;if(Ft>xt*.4&&Tt/Math.abs(_t)>.25||Ft>.9*xt&&xt>a){Q=Math.min(K,o-a),Kt=!0;break}}for(let K=Math.min(s-1,o+a);K<T-l;K++)if(ot[K]*_t<0&&Math.abs(ot[K])>le){let Ft=0,Tt=0,xt=0;for(let wt=K;wt<T;wt++)ot[wt]*_t<0&&(Ft++,Tt=Math.max(Tt,Math.abs(ot[wt]))),xt++;if(Ft>xt*.4&&Tt/Math.abs(_t)>.25||Ft>.9*xt&&xt>a){T=Math.max(K,o+a),Kt=!0;break}}if(Kt&&((o-Q<Math.max(1,Math.round(4/e))||T-o<Math.max(1,Math.round(4/e)))&&(ke=!0),j<2)){j++;continue}return{sampled:Z,fftLeft:Q,fftRight:T,leftMean:tt/ct,rightMean:nt/dt,peakSlopeIdx:pt,slopes:ot,clipped:Kt,dodgy:ke}}})(0,s-1);if(!d)return null;const f=2,p=Math.max(d.leftMean,d.rightMean),g=Math.min(d.leftMean,d.rightMean);let y=d.fftLeft,x=d.fftLeft,b=1/0,_=1/0;for(let E=d.fftLeft;E<=d.fftRight;E++){let J=0,Q=0;for(let z=-f;z<=f;z++){const Z=d.sampled[Math.max(0,Math.min(s-1,E+z))];J+=Z,Q++}const T=J/Math.max(1,Q),j=Math.abs(T-g-.1*(p-g)),G=Math.abs(T-g-.9*(p-g));j<b&&(b=j,y=E),G<_&&(_=G,x=E)}const P=Math.max(4,Math.abs(y-x)*e);if(y<x){const E=y;y=x,x=E}const S=Math.max(l,l+2*Math.round(P/Math.max(e,1e-6)));y+=S,x-=S;const w=Math.max(Math.abs(y-o),Math.abs(x-o),Math.max(a,Math.round(4/Math.max(e,1e-6)))),C=1.85,A=.5,k=new Array(s).fill(0),F=new Array(s).fill(0);for(let E=0;E<n.length;E++){const J=Math.max(0,Math.min(s-1,Math.trunc((n[E]-i)/e)));let Q=5;Math.abs(J-o)>A*w&&(Q=Math.abs(J-o)>2*A*w?12:7);const T=Math.max(d.fftLeft,J-Q),j=Math.min(d.fftRight-1,J+Q);if(j<o-C*w||T>o+C*w){for(let G=T;G<=j;G++)k[G]+=t[E],F[G]+=1;continue}for(let G=T;G<=j;G++){let z=1;if(Math.abs(G-o)<C*w){const Z=c(G);if(Math.abs(G-o)<w*A)z=mn(n[E]-Z,pn,e,1);else{const tt=(Math.abs(G-o)/Math.max(1e-6,w)-A)/Math.max(1e-6,C-A),nt=1*(1-tt)+.01*tt;z=mn(n[E]-Z,pn,e,nt)}}!(z>0)||!Number.isFinite(z)||(k[G]+=t[E]*z,F[G]+=z)}}const M=new Array(s).fill(0);let v=-1,I=-1;for(let E=0;E<s;E++)F[E]>0?(M[E]=k[E]/F[E],v<0&&(v=E),I=E):M[E]=Zt;if(v<0||I<0)return null;const R=Math.max(1,Math.round(3/Math.max(e,1e-6)));let U=M[v],Y=1;for(let E=v+1;E<o&&Y<R;E++)M[E]!==Zt&&(U+=M[E],Y++);U/=Y;let B=M[I],L=1;for(let E=I-1;E>o&&L<R;E--)M[E]!==Zt&&(B+=M[E],L++);B/=L;for(let E=v-1;E>=0;E--)M[E]=U;for(let E=I+1;E<s;E++)M[E]=B;for(let E=v+1;E<I;E++){if(M[E]!==Zt)continue;let J=E-1;for(;J>=0&&M[J]===Zt;)J--;let Q=E+1;for(;Q<s&&M[Q]===Zt;)Q++;if(J>=0&&Q<s){const T=(E-J)/Math.max(1,Q-J);M[E]=M[J]*(1-T)+M[Q]*T}else J>=0?M[E]=M[J]:Q<s&&(M[E]=M[Q])}const D=[...U<=B?M:[...M].reverse()],X=U<=B?D:D.reverse(),et=Math.max(Math.round(o-C*w),d.fftLeft+2),V=Math.min(Math.round(o+C*w),d.fftRight-3),$=Math.max(1,Math.round(2/Math.max(e,1e-6))),H=ys(X,d.fftLeft,d.fftRight,et,V,$),it=new Array(s).fill(0);let at=H[Math.max(0,Math.min(s-1,d.fftLeft))]??H[0]??0;for(let E=d.fftLeft;E<=d.fftRight;E++){const J=H[E]??at,Q=H[Math.min(s-1,E+1)]??J;it[E]=Q-at,at=J}return{esf:H,lsfFull:it}}function bs(n,t,e=Lt){if(n.length===0||t.length!==n.length)return null;const i=se,r=i/2,s=ui,o=2*Et,a=Math.max(1,Math.round(.5*Et)),l=5,u=Math.max(0,Math.round(r-e*Et)),h=Math.min(i-1,Math.round(r+e*Et));if(h-u<32)return null;let c=new Array(i).fill(Zt),m=0,d=0,f=0,p=0,g=u,y=h,x=u,b=h,_=0;for(;;){const T=new Array(i).fill(0),j=new Array(i).fill(0);c=new Array(i).fill(Zt),m=0,d=0,f=0,p=0;let G=-1,z=-1;for(let W=0;W<n.length;W++){const ut=Math.trunc(n[W]*Et+r),ot=Math.max(x,ut-l),pt=Math.min(b-1,ut+l);for(let mt=ot;mt<=pt;mt++){const te=(mt-r)*s,Dt=Math.max(0,1-Math.abs((n[W]-te)*1.75));Dt>0&&(j[mt]+=t[W]*Dt,T[mt]+=Dt)}}const Z=Math.max(r-i/8,x+2*Et),tt=Math.min(r+i/8,b-2*Et);for(let W=Math.max(0,x-1);W<=Math.min(i-1,b+1);W++)T[W]>0&&(c[W]=j[W]/T[W],W<Z&&(m+=c[W],f++),W>tt&&(d+=c[W],p++),G<0&&(G=W),z=W);if(G<0||z<0||f<=0||p<=0)return null;for(let W=G-1;W>=0;W--)c[W]=c[G];for(let W=z+1;W<i;W++)c[W]=c[z];const nt=new Array(i).fill(0);let ct=r;const dt=2*Et;for(let W=dt+1;W<i-1-dt;W++){let ut=0,ot=0;for(let pt=-dt;pt<=dt;pt++)ot+=c[W+pt]*pt,ut+=pt*pt;nt[W]=ut>0?ot/ut:0,Math.abs(nt[W])>Math.abs(nt[ct]??0)&&W>x+dt&&W<b-dt-1&&(ct=W)}if(Math.abs(ct-r)>2*Et&&Math.abs(ct-r)<12*Et)return null;let rt=0;for(let W=Math.max(0,r-dt);W<=Math.min(i-1,r+dt);W++)Math.abs(nt[W])>Math.abs(rt)&&(rt=nt[W]);if(!Number.isFinite(rt)||Math.abs(rt)<=1e-9)return null;const Bt=Math.abs(rt*.001);g=x,y=b;let Vt=!1;for(let W=r-o;W>=x+a;W--)if(nt[W]*rt<0&&Math.abs(nt[W])>Bt){let ut=0,ot=0,pt=0;for(let mt=W;mt>=x;mt--)nt[mt]*rt<0&&(ut++,ot=Math.max(ot,Math.abs(nt[mt]))),pt++;if(ut>pt*.4&&ot/Math.abs(rt)>.25||ut>.9*pt&&pt>o){g=Math.min(W,r-o),Vt=!0;break}}for(let W=r+o;W<b-a;W++)if(nt[W]*rt<0&&Math.abs(nt[W])>Bt){let ut=0,ot=0,pt=0;for(let mt=W;mt<b;mt++)nt[mt]*rt<0&&(ut++,ot=Math.max(ot,Math.abs(nt[mt]))),pt++;if(ut>pt*.4&&ot/Math.abs(rt)>.25||ut>.9*pt&&pt>o){y=Math.max(W,r+o),Vt=!0;break}}if(Vt&&_<2){x=g,b=y,_++;continue}break}const P=Math.max(m/f,d/p),S=Math.min(m/f,d/p);let w=g,C=g,A=1/0,k=1/0;for(let T=g;T<=y;T++){const j=(c[Math.max(0,T-2)]+c[Math.max(0,T-1)]+c[T]+c[Math.min(i-1,T+1)]+c[Math.min(i-1,T+2)])/5,G=Math.abs(j-S-.1*(P-S)),z=Math.abs(j-S-.9*(P-S));G<A&&(A=G,w=T),z<k&&(k=z,C=T)}if(w<C){const T=w;w=C,C=T}const F=Math.max(4,Math.abs(w-C)*s),M=Math.max(a,a+2*Math.trunc(F/Math.max(s,1e-6)));w+=M,C-=M;const v=Math.max(Math.abs(w-r),Math.abs(C-r),Math.max(o,Math.trunc(4/Math.max(s,1e-6)))),I=new Array(i).fill(0),R=new Array(i).fill(0),U=1.85,Y=.5;for(let T=0;T<n.length;T++){const j=Math.trunc(n[T]*Et+r);let G=5;Math.abs(j-r)>Y*v&&(G=Math.abs(j-r)>2*Y*v?12:7);const z=Math.max(g,j-G),Z=Math.min(y-1,j+G);if(Z<r-U*v||z>r+U*v){for(let tt=z;tt<=Z;tt++)I[tt]+=t[T],R[tt]+=1;continue}for(let tt=z;tt<=Z;tt++){let nt=1;if(Math.abs(tt-r)<U*v){const ct=(tt-r)*s;if(Math.abs(tt-r)<v*Y)nt=mn(n[T]-ct,pn,s,1);else{const dt=(Math.abs(tt-r)/Math.max(1e-6,v)-Y)/Math.max(1e-6,U-Y),rt=1*(1-dt)+.01*dt;nt=mn(n[T]-ct,pn,s,rt)}}!(nt>0)||!Number.isFinite(nt)||(I[tt]+=t[T]*nt,R[tt]+=nt)}}const B=new Array(i).fill(0);let L=-1,O=-1;for(let T=Math.max(0,g-1);T<=Math.min(i-1,y+1);T++)R[T]>0?(B[T]=I[T]/R[T],L<0&&(L=T),O=T):B[T]=Zt;if(L<0||O<0)return null;const D=3*Et;let X=B[L],et=1;for(let T=L+1;T<r&&et<D;T++)B[T]!==Zt&&(X+=B[T],et++);X/=Math.max(1,et);let V=B[O],$=1;for(let T=O-1;T>r&&$<D;T--)B[T]!==Zt&&(V+=B[T],$++);V/=Math.max(1,$);for(let T=L-1;T>=0;T--)B[T]=X;for(let T=O+1;T<i;T++)B[T]=V;const H=Math.max(Math.trunc(r-U*v),g+2),it=Math.min(Math.trunc(r+U*v),y-3),at=Math.max(1,Math.trunc(2/Math.max(s,1e-6))),E=ys(B,g,y,H,it,at),J=new Array(i).fill(0);let Q=E[Math.max(0,Math.min(i-1,g))]??E[0]??0;for(let T=g;T<=y;T++){const j=E[T]??Q;J[T]=(E[Math.min(i-1,T+1)]??j)-Q,Q=j}return{esf:E,lsfFull:J}}function yo(n){const t=new Array(n.length).fill(0);if(n.length===0)return t;t[0]=n[0];for(let e=1;e<n.length;e++){const i=ir(n[e]-n[e-1]);t[e]=t[e-1]+i}return t}function ir(n){if(!Number.isFinite(n))return n;let t=(n+Math.PI)%(2*Math.PI);return t<0&&(t+=2*Math.PI),t-Math.PI}function xo(n,t,e=0){if(n.length===0)return[];const i=Number.isFinite(e)?e:0,r=n.map((o,a)=>{const l=t[a]??0,u=-2*Math.PI*i*l;return ir(o-u)});return yo(r).map((o,a)=>{const l=t[a]??0,u=-2*Math.PI*i*l;return o+u})}function bo(n,t,e,i=.05,r=Number.POSITIVE_INFINITY){const s=Math.min(n.length,t.length);if(s<2)return null;let o=0;if(e)for(let p=1;p<s;p++){const g=e[p];Number.isFinite(g)&&g>o&&(o=g)}const a=e&&o>0?Math.max(1e-6,o*i):0;let l=0,u=0,h=0,c=0,m=0,d=0;for(let p=1;p<s;p++){const g=n[p],y=t[p];if(!Number.isFinite(g)||!Number.isFinite(y)||Math.abs(g)<=1e-12||g>r)continue;const x=e?e[p]:1;if(!Number.isFinite(x)||x<=a)continue;const b=e?x*x:1;l+=b,u+=b*g,h+=b*y,c+=b*g*g,m+=b*g*y,d++}if(d<4||l<=0)return null;const f=l*c-u*u;return Math.abs(f)<=1e-12?null:{slope:(l*m-u*h)/f,intercept:(h*c-u*m)/f,used:d,threshold:a}}function _o(n,t,e=Number.POSITIVE_INFINITY){const i=[],r=[],s=Math.min(n.length,t.length);for(let o=1;o<s;o++)Number.isFinite(n[o])&&Number.isFinite(t[o])&&Math.abs(n[o])>1e-12&&n[o]<=e&&(i.push(n[o]),r.push(t[o]));return i.length<2?{slope:0,intercept:Number.isFinite(t[0])?t[0]:0,used:i.length}:{...nr(i,r),used:i.length}}function Mo(n,t,e,i=.05,r=Number.POSITIVE_INFINITY,s=0){const o=Math.min(n.length,t.length);if(o<4)return null;let a=0;if(e)for(let b=1;b<o;b++){const _=e[b];Number.isFinite(_)&&_>a&&(a=_)}const l=e&&a>0?Math.max(1e-6,a*i):0,u=[];for(let b=1;b<o;b++){const _=n[b],P=t[b];if(!Number.isFinite(_)||!Number.isFinite(P)||Math.abs(_)<=1e-12||_>r)continue;const S=e?e[b]:1;!Number.isFinite(S)||S<=l||u.push({freq:_,phase:P,weight:e?S*S:1})}if(u.length<4)return null;const h=b=>{let _=0,P=0;for(const k of u){const F=k.phase+2*Math.PI*b*k.freq;_+=k.weight*Math.sin(F),P+=k.weight*Math.cos(F)}const S=Math.atan2(_,P);let w=0,C=0;const A=.65;for(const k of u){const F=k.phase+2*Math.PI*b*k.freq,M=Math.abs(ir(F-S)),v=M<=A?M*M:A*(2*M-A);w+=k.weight*v,C+=k.weight}return{score:C>0?w/C:Number.POSITIVE_INFINITY,intercept:S}},c=Number.isFinite(s)?s:0,m=Math.max(2,Math.min(8,Math.abs(c)>1e-6?4:2)),d=.02;let f=c,p=h(f);for(let b=c-m;b<=c+m+d*.5;b+=d){const _=h(b);_.score<p.score&&(p=_,f=b)}let g=f-d*2,y=f+d*2;for(let b=0;b<32;b++){const _=g+(y-g)/3,P=y-(y-g)/3,S=h(_).score,w=h(P).score;S<w?y=P:g=_}const x=(g+y)*.5;return p=h(x),{slope:-2*Math.PI*x,intercept:p.intercept,used:u.length,threshold:l}}function wo(n,t){if(Number.isFinite(t)&&t>0)return t;let e=0;for(const i of n)Number.isFinite(i)&&i>e&&(e=i);return e>0?e:Number.POSITIVE_INFINITY}function _s(n,t,e,i=Number.POSITIVE_INFINITY,r=0){const s=xo(n,t,r),o=wo(t,i),a=Mo(t,n,e,.05,o,r),l=a?null:bo(t,s,e,.05,o),u=a||l?null:_o(t,s,o),h=(a==null?void 0:a.slope)??(l==null?void 0:l.slope)??(u==null?void 0:u.slope)??0,c=(a==null?void 0:a.intercept)??(l==null?void 0:l.intercept)??(u==null?void 0:u.intercept)??0,m=s.map((x,b)=>x-(h*(t[b]??0)+c)),d=Number.isFinite(m[0])?m[0]:0,f=t.map(x=>h*x+c+d),p=m.map(x=>x-d),g=c+d,y=Number.isFinite(h)?-h/(2*Math.PI):null;return{raw:[...n],unwrapped:s,linear:f,residual:p,fit:{groupDelayPx:y===null?null:y-r,absoluteGroupDelayPx:y,referenceDelayPx:r,slopeRadPerCycle:Number.isFinite(h)?h:null,interceptRad:Number.isFinite(g)?g:null,fitPointCount:(a==null?void 0:a.used)??(l==null?void 0:l.used)??(u==null?void 0:u.used)??0,fitWeightThreshold:(a==null?void 0:a.threshold)??(l==null?void 0:l.threshold)??0,fitDomain:"cycles-per-pixel",fitMaxFreqCyclesPerPixel:o}}}function Ms(n,t,e){const i=[],r=[],s=[],o=[];for(let a=0;a<e.length;a++){const l=e[a];i.push(Ze(l,t,n.raw)),r.push(Ze(l,t,n.unwrapped)),s.push(Ze(l,t,n.linear)),o.push(Ze(l,t,n.residual))}return{ptfRaw:i,ptfUnwrapped:r,ptfLinear:s,ptfResidual:o}}function ws(n,t){const e=n.map((a,l)=>({dist:a,value:t[l]})).filter(a=>Number.isFinite(a.dist)&&Number.isFinite(a.value)).sort((a,l)=>a.dist-l.dist);if(e.length===0)return{dists:[],vals:[]};const i=Math.max(1,Math.min(16,Math.floor(e.length*.1)));let r=0,s=0;for(let a=0;a<i;a++)r+=e[a].value,s+=e[e.length-1-a].value;r/=i,s/=i;const o=r<=s?e:e.map(a=>({dist:-a.dist,value:a.value})).sort((a,l)=>a.dist-l.dist);return{dists:o.map(a=>a.dist),vals:o.map(a=>a.value)}}function So(n,t,e="RGGB"){const i=e.toUpperCase(),r=t%2,s=n%2;return i==="RGGB"||i==="BGGR"?(r+s)%2!==0:i==="GBRG"||i==="GRBG"?(r+s)%2===0:(r+s)%2!==0}function vo(n,t){return n+t&1?2:1}function Or(n,t){return(t&1)<<1|n&1}function xn(n,t,e){return n===void 0?0:typeof n=="number"?Number.isFinite(n)?n:0:Number.isFinite(n[Or(t,e)])?n[Or(t,e)]:0}function vt(n,t,e,i){return i!==void 0&&i!=="default"?vo(n,t)===i:So(n,t,e)}function ei(n){return n.length===0?0:n.reduce((t,e)=>t+e,0)/n.length}function Co(n,t,e){const i=(e%t+t)%t,r=Math.floor(i),s=(r+1)%t,o=i-r,a=r<n.length?n[r]:0,l=s<n.length?n[s]:0;return a*(1-o)+l*o}function Ss(n,t){const e=n.length;if(e===0)return 0;if(t<=0)return n[0];if(t>=e-1)return n[e-1];const i=Math.floor(t),r=Math.min(e-1,i+1),s=t-i;return n[i]*(1-s)+n[r]*s}function Po(n,t,e){const i=n.length,r=new Array(i).fill(0);for(let s=0;s<i;s++)r[s]=Ss(n,s-e+t);return r}function ni(n,t,e){const i=n.length;if(i===0)return{peakPos:0,peakIdx:0,peakVal:0};const r=Math.max(0,Math.floor(t-e)),s=Math.min(i-1,Math.ceil(t+e));let o=Math.max(0,Math.min(i-1,Math.round(t))),a=-1/0;for(let u=r;u<=s;u++){const h=Math.abs(n[u]);h>a&&(a=h,o=u)}Number.isFinite(a)||(a=Math.abs(n[o]??0));let l=o;if(o>0&&o<i-1){const u=n[o]>=0?1:-1,h=u*n[o-1],c=u*n[o],m=u*n[o+1],d=h-2*c+m;if(Number.isFinite(d)&&Math.abs(d)>1e-9){const f=.5*(h-m)/d;Number.isFinite(f)&&Math.abs(f)<=1&&(l=o+f)}}return{peakPos:l,peakIdx:o,peakVal:Math.abs(Ss(n,l))}}function ko(n,t,e,i,r,s,o,a,l){const u=Math.floor(i.x),h=Math.floor(i.y),c=Math.floor(i.w),m=Math.floor(i.h),d=[],f=(p,g)=>{if(p<0||g<0||p>=t||g>=e)return null;const y=r+p,x=s+g;return Math.max(0,n[g*t+p]-xn(l,y,x))};for(let p=0;p<m;p++){const g=[],y=h+p;for(let x=0;x<c;x++){const b=u+x,_=r+b,P=s+y;if(vt(_,P,o,a)){g.push(f(b,y)??0);continue}const S=[],w=f(b-1,y),C=f(b+1,y),A=f(b,y-1),k=f(b,y+1);if(w!==null&&vt(_-1,P,o,a)&&S.push(w),C!==null&&vt(_+1,P,o,a)&&S.push(C),A!==null&&vt(_,P-1,o,a)&&S.push(A),k!==null&&vt(_,P+1,o,a)&&S.push(k),S.length===0){const F=[],M=f(b-1,y-1),v=f(b+1,y-1),I=f(b-1,y+1),R=f(b+1,y+1);M!==null&&vt(_-1,P-1,o,a)&&F.push(M),v!==null&&vt(_+1,P-1,o,a)&&F.push(v),I!==null&&vt(_-1,P+1,o,a)&&F.push(I),R!==null&&vt(_+1,P+1,o,a)&&F.push(R),g.push(ei(F));continue}g.push(ei(S))}d.push(g)}return d}function vs(n,t,e,i,r,s,o,a,l){const u=Math.floor(i.x),h=Math.floor(i.y),c=Math.floor(i.w),m=Math.floor(i.h),d=[];for(let f=0;f<m;f++){const p=h+f;for(let g=0;g<c;g++){const y=u+g,x=r+y,b=s+p;vt(x,b,o,a)&&d.push({x,y:b,value:Math.max(0,n[p*t+y]-xn(l,x,b))})}}return d}function Cs(n,t,e,i,r){const s=Math.floor(i.x),o=Math.floor(i.y),a=Math.floor(i.w),l=Math.floor(i.h),u=(r==null?void 0:r.globalX)??0,h=(r==null?void 0:r.globalY)??0,c=!!(r!=null&&r.isThreePlane)&&n.length>=t*e*3,m=r==null?void 0:r.threePlaneChannel,d=[];for(let f=0;f<l;f++){const p=o+f,g=p*t;for(let y=0;y<a;y++){const x=s+y;let b=0;if(!c)b=Math.max(0,n[g+x]-xn(r==null?void 0:r.blackLevel,u+x,h+p));else{const _=(g+x)*3;if(m!==void 0)b=n[_+m];else{const P=n[_],S=n[_+1],w=n[_+2];b=.2126*P+.7152*S+.0722*w}}d.push({x:u+x,y:h+p,value:b})}}return d}function Fo(n){var s;const t=n.length,e=((s=n[0])==null?void 0:s.length)??0;let i=0,r=0;for(let o=1;o<t-1;o++)for(let a=1;a<e-1;a++)i+=Math.abs(n[o][a+1]-n[o][a-1]),r+=Math.abs(n[o+1][a]-n[o-1][a]);return{gx:i,gy:r}}function Ps(n,t,e,i,r,s,o){var m;const a=n.length,l=((m=n[0])==null?void 0:m.length)??0,u=(d,f,p,g,y,x,b)=>{const _=b?a:l,P=Math.max(0,f-3),S=Math.min(_,f+4);let w=0,C=0;for(let k=P;k<S;k++)w+=d[k],C+=k*d[k];if(w<=0)return null;const A=C/w;return b?{x:t+g*y,y:p+A*x,weight:w}:{x:p+A*x,y:e+g*y,weight:w}},h=(d,f,p,g,y,x,b)=>{const _=Math.max(3,Math.min(Math.max(3,Math.floor(p/3)),Math.max(4,Math.round(p*.12)))),P=d.map(v=>{let I=-1/0,R=-1;for(let U=0;U<v.length;U++)v[U]>I&&(I=v[U],R=U);return{peakValue:I,peakIndex:R}}),S=(f-1)*.5,w=P.map((v,I)=>({...v,index:I})).filter(v=>v.peakValue>1&&v.peakIndex>=0).sort((v,I)=>{const R=I.peakValue-v.peakValue;return Math.abs(R)>1e-6?R:Math.abs(v.index-S)-Math.abs(I.index-S)});if(w.length===0)return[];const C=w[0],A=new Array(f).fill(null),k=u(d[C.index],C.peakIndex,g,C.index,y,x,b);if(!k)return[];A[C.index]=k;const F=(v,I)=>{const R=d[v],U=P[v];if(!(U.peakValue>1)||U.peakIndex<0)return null;const Y=Math.max(0,Math.floor(I-_)),B=Math.min(R.length,Math.ceil(I+_+1));let L=-1/0,O=-1;for(let V=Y;V<B;V++)R[V]>L&&(L=R[V],O=V);if(O<0||!(L>1))return null;const D=Math.max(1e-6,U.peakValue),X=Math.abs(O-I)<=_,et=L>=D*.25;return!X||!et?null:u(R,O,g,v,y,x,b)};let M=k?b?(k.y-g)/x:(k.x-g)/x:C.peakIndex;for(let v=C.index+1;v<f;v++){const I=F(v,M);I&&(A[v]=I,M=b?(I.y-g)/x:(I.x-g)/x)}M=k?b?(k.y-g)/x:(k.x-g)/x:C.peakIndex;for(let v=C.index-1;v>=0;v--){const I=F(v,M);I&&(A[v]=I,M=b?(I.y-g)/x:(I.x-g)/x)}return A.filter(v=>!!v)};if(s){const d=n.map(f=>f.map((p,g)=>g===0?0:Math.abs(p-f[g-1])));return h(d,a,l,t,r,i,!1)}const c=[];for(let d=0;d<l;d++){const f=new Array(a).fill(0);for(let p=1;p<a;p++)f[p]=Math.abs(n[p][d]-n[p-1][d]);c.push(f)}return h(c,l,a,e,i,r,!0)}function ve(n){if(n.length<2)return null;let t=0,e=0,i=0;for(const c of n)t+=c.weight,e+=c.x*c.weight,i+=c.y*c.weight;if(t<=0)return null;e/=t,i/=t;let r=0,s=0,o=0;for(const c of n){const m=c.x-e,d=c.y-i;r+=c.weight*m*m,s+=c.weight*d*d,o+=c.weight*m*d}r/=t,s/=t,o/=t;const a=.5*Math.atan2(2*o,r-s);let l=Math.cos(a),u=Math.sin(a);const h=Math.hypot(l,u);return!Number.isFinite(h)||h<=1e-9?null:(l/=h,u/=h,(l<0||Math.abs(l)<=1e-9&&u<0)&&(l=-l,u=-u),{pointX:e,pointY:i,dirX:l,dirY:u,orientation:Math.abs(l)>=Math.abs(u)?1:2})}function Ao(n,t){if(n.length!==4||t.length!==4||n.some(i=>i.length!==4))return null;const e=n.map((i,r)=>[...i,t[r]]);for(let i=0;i<4;i++){let r=i,s=Math.abs(e[i][i]);for(let a=i+1;a<4;a++){const l=Math.abs(e[a][i]);l>s&&(s=l,r=a)}if(!(s>1e-12))return null;if(r!==i){const a=e[i];e[i]=e[r],e[r]=a}const o=e[i][i];for(let a=i;a<=4;a++)e[i][a]/=o;for(let a=0;a<4;a++){if(a===i)continue;const l=e[a][i];if(!(Math.abs(l)<=1e-12))for(let u=i;u<=4;u++)e[a][u]-=l*e[i][u]}}return[e[0][4],e[1][4],e[2][4],e[3][4]]}function To(n){if(n.length<4)return 0;const t=[...n].sort((m,d)=>m.x-d.x),e=t[0].x,r=t[t.length-1].x-e;if(!(r>1e-6))return 0;const s=16,o=[];for(let m=0;m<s;m++){const d=Math.max(0,Math.floor((m-1.5)*t.length/s)),f=Math.min(t.length-1,Math.floor((m+2.5)*t.length/s));if(f<d)continue;let p=0,g=0,y=0;for(let x=d;x<=f;x++)p+=t[x].x,g+=t[x].y,y++;y>0&&o.push({x:p/y,y:g/y})}if(o.length<4)return 0;const a=[.05952381,0,-.03571429,-.04761905,-.03571429,0,.05952381],l=new Array(o.length).fill(0),u=3;for(let m=0;m<o.length;m++){let d=0;for(let f=-u;f<=u;f++){const p=m+f;p<0||p>=o.length||(d+=a[f+u]*o[p].y)}l[m]=d}let h=0,c=1/0;for(let m=0;m<l.length-1;m++){const d=l[m],f=l[m+1];if(d===0){const b=Math.abs(o[m].x);b<c&&(c=b,h=o[m].x);continue}if(d*f>=0)continue;const p=f-d,g=Math.abs(p)>1e-12?-d/p:.5,y=o[m].x+(o[m+1].x-o[m].x)*g,x=Math.abs(y);x<c&&(c=x,h=y)}return!Number.isFinite(h)||h<e+.3*r||h>e+.7*r?0:h}function Io(n){if(n.length<8)return null;const t=[...n].filter(f=>Number.isFinite(f.x)&&Number.isFinite(f.y)&&Number.isFinite(f.weight)).sort((f,p)=>f.x-p.x);if(t.length<8)return null;const i=[To(t),0,.5*(t[Math.floor((t.length-1)*.5)].x+t[Math.ceil((t.length-1)*.5)].x)];let r=null;for(const f of i){if(!Number.isFinite(f))continue;let p=0,g=0;for(const y of t)y.x<=f?p++:g++;if(p>=4&&g>=4){r=f;break}}if(r===null)return null;const s=Array.from({length:4},()=>new Array(4).fill(0)),o=new Array(4).fill(0);for(const f of t){const p=f.x,g=f.y,y=Math.max(1e-6,f.weight),x=p<=r?[p*p,p,1,0]:[2*r*p-r*r,p,1,(p-r)*(p-r)];for(let b=0;b<4;b++){o[b]+=y*x[b]*g;for(let _=0;_<4;_++)s[b][_]+=y*x[b]*x[_]}}const a=Ao(s,o);if(!a)return null;const[l,u,h,c]=a,m=u+2*(l-c)*r,d=h+(c-l)*r*r;return[l,u,h,c,m,d].every(f=>Number.isFinite(f))?{splitX:r,left:[l,u,h],right:[c,m,d]}:null}function No(n,t,e){const[i,r,s]=e;if(Math.abs(i)<=1e-12){const b=1+r*r;return b>1e-12?[(n-r*(s-t))/b]:[n]}const o=2*i*i,a=3*i*r,l=1+2*i*s-2*i*t+r*r,u=r*s-t*r-n;if(Math.abs(o)<=1e-12)return[n];const h=a/o,c=l/o,m=u/o,d=(h*h-3*c)/9,f=(2*h*h*h-9*h*c+27*m)/54,p=f*f-d*d*d;if(p<0&&d>0){const b=Math.acos(Math.max(-1,Math.min(1,f/Math.sqrt(d*d*d)))),_=-2*Math.sqrt(d);return[_*Math.cos(b/3)-h/3,_*Math.cos((b+2*Math.PI)/3)-h/3,_*Math.cos((b-2*Math.PI)/3)-h/3]}const g=Math.sqrt(Math.max(0,p)),y=-Math.sign(f||1)*Math.cbrt(Math.abs(f)+g),x=Math.abs(y)<=1e-12?0:d/y;return[y+x-h/3]}function ks(n,t){if(n.length<8)return null;const e=-t.dirY,i=t.dirX,r=n.map(o=>({x:(o.x-t.pointX)*t.dirX+(o.y-t.pointY)*t.dirY,y:(o.x-t.pointX)*e+(o.y-t.pointY)*i,weight:o.weight})),s=Io(r);return s?{...t,normalX:e,normalY:i,splitX:s.splitX,left:s.left,right:s.right}:null}function Ro(n,t){const e=n.x-t.pointX,i=n.y-t.pointY,r=e*t.dirX+i*t.dirY,s=e*t.normalX+i*t.normalY,o=r<t.splitX?t.left:t.right,a=No(r,s,o);let l=s,u=Number.POSITIVE_INFINITY;for(const h of a){if(!Number.isFinite(h))continue;const c=o[0]*h*h+o[1]*h+o[2],m=r-h,d=s-c,f=Math.hypot(m,d);Number.isFinite(f)&&f<u&&(u=f,l=(d>=0?1:-1)*f)}return Number.isFinite(u)?l:s}function Fs(n,t,e,i,r,s,o,a=Lt){if(!t||t.length<8||n.length===0)return null;const l=t.filter(v=>Number.isFinite(v.x)&&Number.isFinite(v.y)).map(v=>({x:v.x,y:v.y,weight:1}));if(l.length<8)return null;const u=ve(l);if(!u)return null;const h=e.p2.x-e.p1.x,c=e.p2.y-e.p1.y,m=Math.hypot(h,c);if(!Number.isFinite(m)||m<=1e-6)return null;let d=u.dirX,f=u.dirY;d*h+f*c<0&&(d=-d,f=-f);const p={...u,dirX:d,dirY:f},g=ks(l,p),y=h/m,x=c/m,b=-x,_=y,P=(e.p1.x+e.p2.x)*.5,S=(e.p1.y+e.p2.y)*.5,w=-p.dirY,C=p.dirX,A=Math.abs(y)>=Math.abs(x)?1:2,k=[],F=[];for(const v of n){const I=v.x-P,R=v.y-S,U=I*y+R*x;if(Math.abs(U)>i)continue;const Y=I*b+R*_;if(Math.abs(Y)>r)continue;const B=g?Ro(v,g):(v.x-p.pointX)*w+(v.y-p.pointY)*C;Number.isFinite(B)&&(k.push(B),F.push(v.value))}if(k.length<8)return null;const M=o?s!=null&&s.forceLegacyModel?An(k,F,A,a,r*2):bn(k,F,A,a):hi(k,F,Math.max(2,r*2),s==null?void 0:s.manualBinSize,A,s==null?void 0:s.preferAutoPerEdgeBin);return M?(M.quadraticProjectionUsed=!!g,M):null}function Gr(n){if(n.length<2)return null;const t=n.filter(e=>Number.isFinite(e.x)&&Number.isFinite(e.y)).map(e=>({x:e.x,y:e.y,weight:1}));return t.length<2?null:Cn(t,ve(t))}function Cn(n,t){if(!t||n.length<2)return null;let e=1/0,i=-1/0;for(const s of n){const o=(s.x-t.pointX)*t.dirX+(s.y-t.pointY)*t.dirY;e=Math.min(e,o),i=Math.max(i,o)}if(!Number.isFinite(e)||!Number.isFinite(i))return null;const r=Math.max(.5,(i-e)*.03);return{p1:{x:t.pointX+t.dirX*(e-r),y:t.pointY+t.dirY*(e-r)},p2:{x:t.pointX+t.dirX*(i+r),y:t.pointY+t.dirY*(i+r)}}}function Lo(n,t,e,i,r,s,o,a){return[Cn(n,t),Cn(e,i),Cn(r,s),Cn(o,a)]}function Eo(n,t,e){if(!n||n.length<8)return;const i=n.filter(x=>Number.isFinite(x.x)&&Number.isFinite(x.y)).map(x=>({x:x.x,y:x.y,weight:1}));if(i.length<8)return;const r=ve(i);if(!r)return;const s=t.p2.x-t.p1.x,o=t.p2.y-t.p1.y,a=Math.hypot(s,o);if(!Number.isFinite(a)||a<=1e-6)return;let l=r.dirX,u=r.dirY;l*s+u*o<0&&(l=-l,u=-u);const h=ks(i,{...r,dirX:l,dirY:u});if(!h)return;const c=n.map(x=>(x.x-h.pointX)*h.dirX+(x.y-h.pointY)*h.dirY).filter(x=>Number.isFinite(x)),m=(t.p1.x-h.pointX)*h.dirX+(t.p1.y-h.pointY)*h.dirY,d=(t.p2.x-h.pointX)*h.dirX+(t.p2.y-h.pointY)*h.dirY;if(Number.isFinite(m)&&c.push(m),Number.isFinite(d)&&c.push(d),c.length<2)return;const f=Math.min(...c),p=Math.max(...c);if(!Number.isFinite(f)||!Number.isFinite(p)||p-f<=1e-6)return;const g=Math.max(21,e),y=[];for(let x=0;x<g;x++){const b=g===1?.5:x/(g-1),_=f+(p-f)*b,P=_<h.splitX?h.left:h.right,S=P[0]*_*_+P[1]*_+P[2];y.push({x:h.pointX+_*h.dirX+S*h.normalX,y:h.pointY+_*h.dirY+S*h.normalY})}return y}function hi(n,t,e,i,r,s=!1,o=!1,a=!1){if(n.length===0||t.length!==n.length)return null;const l=ws(n,t),u=l.dists,h=l.vals;if(u.length===0)return null;const c=()=>{const w=e/2;let C=0;for(const k of n)Math.abs(k)<=w&&C++;if(C<=0)return .125;const A=40*w/C;return Math.max(.01,Math.min(.125,A))},m=(w,C,A,k,F)=>{if(!(k>0)||!(F>0)||!(A>C))return!1;const M=Math.floor((A-C)/k);if(M<2)return!1;const v=Math.max(C,-F),I=Math.min(A,F);if(!(I>v))return!1;const R=Math.max(0,Math.floor((v-C)/k)),U=Math.min(M,Math.ceil((I-C)/k));if(U<=R)return!1;const Y=new Array(U-R).fill(0),B=C+R*k,L=C+U*k;for(let O=0;O<w.length;O++){const D=w[O];if(D<B)continue;if(D>=L)break;const X=Math.floor((D-C)/k);X>=R&&X<U&&Y[X-R]++}return Y.every(O=>O>0)},d=()=>{const w=u[0],C=u[u.length-1],A=Math.max(0,e*.25),k=.125,F=.5,M=.001,v=Math.round((F-k)/M);for(let I=0;I<=v;I++){const R=Number((k+I*M).toFixed(3));if(m(u,w,C,R,A))return R}return F};let f=.125;i&&i>0?f=Math.max(.01,Math.min(.5,i)):s?f=d():f=c();const p=u[0],g=u[u.length-1],y=Math.floor((g-p)/f);if(y<2)return null;const x=()=>{const w=new Array(y).fill(0),C=new Array(y).fill(0);for(let k=0;k<u.length;k++){const F=(u[k]-p)/f;if(Number.isFinite(F))if(o){const M=Math.floor(F),v=F-M,I=1-v,R=v;M>=0&&M<y&&(w[M]+=h[k]*I,C[M]+=I);const U=M+1;U>=0&&U<y&&(w[U]+=h[k]*R,C[U]+=R)}else{const M=Math.floor(F);M>=0&&M<y&&(w[M]+=h[k],C[M]++)}}let A=h[0];for(let k=0;k<y;k++)C[k]>0?(w[k]/=C[k],A=w[k]):w[k]=A;return w},b=a?null:go(u,h,f),_=(b==null?void 0:b.esf)??x(),P=(b==null?void 0:b.lsfFull)??uo(_,f),S=Math.max(0,Math.min(y-1,-p/f-.5));return{esf:_,lsfFull:P,binSize:f,orientation:r,zeroIndex:S,shortSidePx:e,fallbackUsed:a||!b}}function bn(n,t,e,i=Lt){if(n.length===0||t.length!==n.length)return null;const r=ws(n,t),s=r.dists.map((a,l)=>({dist:a,value:r.vals[l]})).filter(a=>Math.abs(a.dist)<i);if(s.length<8)return null;const o=bs(s.map(a=>a.dist),s.map(a=>a.value),i);return o?{esf:o.esf,lsfFull:o.lsfFull,binSize:ui,orientation:e,zeroIndex:se/2,shortSidePx:i*2,fallbackUsed:!1,mtfmapperLike:!0,mtfmapperOrderedDists:s.map(a=>a.dist),mtfmapperOrderedVals:s.map(a=>a.value),mtfmapperEffectiveMaxDot:i}:null}function An(n,t,e,i=Lt,r=i*2){if(n.length===0||t.length!==n.length)return null;const s=[],o=[];for(let a=0;a<n.length;a++){const l=n[a],u=t[a];!Number.isFinite(l)||!Number.isFinite(u)||Math.abs(l)>=i||(s.push(l),o.push(u))}return s.length<8?null:hi(s,o,Math.max(2,r),void 0,e,!0,!0,!0)}function Se(n,t,e,i,r,s,o,a,l=Lt){if(s<=0||o<=0)return null;const h=!(!!(a!=null&&a.isThreePlane)&&n.length>=t*e*3)&&((a==null?void 0:a.greenOnly)??!1),c=(a==null?void 0:a.bayerPattern)||"RGGB",m=r.p2.x-r.p1.x,d=r.p2.y-r.p1.y,f=Math.hypot(m,d);if(!Number.isFinite(f)||f<=1e-6)return null;const p=m/f,g=d/f,y=-g,x=p,b=(r.p1.x+r.p2.x)*.5,_=(r.p1.y+r.p2.y)*.5,P=Math.abs(p)>=Math.abs(g)?1:2,S=h?vs(n,t,e,i,0,0,c,a==null?void 0:a.greenPhase,a==null?void 0:a.blackLevel):Cs(n,t,e,i,{...a,globalX:0,globalY:0});if(S.length===0)return null;if(!(a!=null&&a.disableQuadraticProjection)){const A=Fs(S,a==null?void 0:a.quadraticFitPoints,r,s,o,a,!0,l);if(A)return A}const w=[],C=[];for(const A of S){const k=A.x-b,F=A.y-_,M=k*p+F*g;if(Math.abs(M)>s)continue;const v=k*y+F*x;Math.abs(v)>o||(w.push(v),C.push(A.value))}return w.length<8?null:a!=null&&a.forceLegacyModel?An(w,C,P,l):bn(w,C,P,l)}function Uo(n,t,e=0){const i=[...n.lsfFull];if(i.length<3)return!1;const r=Math.max(n.binSize,1e-6),s=Number.isFinite(n.zeroIndex)?n.zeroIndex:i.length/2,o=Math.max(1,Math.round((n.shortSidePx??0)*.5/r));let{peakPos:a,peakIdx:l,peakVal:u}=ni(i,s,o);const h=u*.2;let c=0,m=i.length-1;for(let p=l;p>=0;p--)if(i[p]<h){c=p;break}for(let p=l;p<i.length;p++)if(i[p]<h){m=p;break}const d=m-c;if(t&&d>0){const p=d*4,g=[],y=[];if(e>0){const x=Math.max(0,l-p-e),b=Math.max(0,l-p);for(let S=x;S<b;S++)g.push(S),y.push(i[S]);const _=Math.min(i.length,l+p),P=Math.min(i.length,l+p+e);for(let S=_;S<P;S++)g.push(S),y.push(i[S])}else{for(let x=0;x<Math.max(0,l-p);x++)g.push(x),y.push(i[x]);for(let x=Math.min(i.length,l+p);x<i.length;x++)g.push(x),y.push(i[x])}if(g.length>2){const{slope:x,intercept:b}=nr(g,y);for(let _=0;_<i.length;_++)i[_]=i[_]-(x*_+b);({peakPos:a}=ni(i,s,o))}}return Math.abs(a-s)*r<=Math.max(1e-6,(n.shortSidePx??0)/6)}function Bo(n){const t=n.length;if(t<3)return!1;let e=0,i=-1/0;for(let o=0;o<t;o++){const a=Math.abs(n[o]);a>i&&(i=a,e=o)}const r=t/3,s=2*t/3;return e>=r&&e<=s}function Mt(n,t,e){const i=Math.max(0,Math.floor(n.x)),r=Math.max(0,Math.floor(n.y)),s=Math.min(t,Math.ceil(n.x+n.w)),o=Math.min(e,Math.ceil(n.y+n.h)),a=s-i,l=o-r;return a<2||l<2?null:{x:i,y:r,w:a,h:l}}function rr(n,t,e,i){const r=[],s=n.x,o=n.y,a=n.x+n.w,l=n.y+n.h,u=n.x+n.w*.5,h=n.y+n.h*.5,c=[{x:s,y:o},{x:a,y:o},{x:a,y:l},{x:s,y:l},{x:u,y:o},{x:a,y:h},{x:u,y:l},{x:s,y:h},{x:u,y:h}];for(const m of c){const d=Pe(m,t);Number.isFinite(d.x)&&Number.isFinite(d.y)&&r.push(d)}return r.length===0?null:Mt(Gt(r,2),e,i)}function Gt(n,t=0){let e=1/0,i=1/0,r=-1/0,s=-1/0;for(const o of n)e=Math.min(e,o.x),i=Math.min(i,o.y),r=Math.max(r,o.x),s=Math.max(s,o.y);return{x:e-t,y:i-t,w:r-e+t*2,h:s-i+t*2}}function Vr(n,t){let e=Math.atan2(t,n)*180/Math.PI;return e<0&&(e+=180),e}function xe(n,t){const e=n.p2.x-n.p1.x,i=n.p2.y-n.p1.y,r=Math.hypot(e,i);if(!Number.isFinite(r)||r<=1e-6)return null;const s=-i/r,o=e/r;return[{x:n.p1.x+s*t,y:n.p1.y+o*t},{x:n.p2.x+s*t,y:n.p2.y+o*t},{x:n.p2.x-s*t,y:n.p2.y-o*t},{x:n.p1.x-s*t,y:n.p1.y-o*t}]}function Do(n,t,e,i,r,s,o,a){if(s<=0||o<=0)return null;const u=!(!!(a!=null&&a.isThreePlane)&&n.length>=t*e*3)&&((a==null?void 0:a.greenOnly)??!1),h=(a==null?void 0:a.bayerPattern)||"RGGB",c=r.p2.x-r.p1.x,m=r.p2.y-r.p1.y,d=Math.hypot(c,m);if(!Number.isFinite(d)||d<=1e-6)return null;const f=c/d,p=m/d,g=-p,y=f,x=(r.p1.x+r.p2.x)*.5,b=(r.p1.y+r.p2.y)*.5,_=Math.abs(f)>=Math.abs(p)?1:2,P=u?vs(n,t,e,i,0,0,h,a==null?void 0:a.greenPhase,a==null?void 0:a.blackLevel):Cs(n,t,e,i,{...a,globalX:0,globalY:0});if(P.length===0)return null;if(!(a!=null&&a.disableQuadraticProjection)){const C=Fs(P,a==null?void 0:a.quadraticFitPoints,r,s,o,a,!1);if(C)return C}const S=[],w=[];for(const C of P){const A=C.x-x,k=C.y-b,F=A*f+k*p;if(Math.abs(F)>s)continue;const M=A*g+k*y;Math.abs(M)>o||(S.push(M),w.push(C.value))}return S.length<8?null:hi(S,w,o*2,a==null?void 0:a.manualBinSize,_,a==null?void 0:a.preferAutoPerEdgeBin)}function Oo(n,t,e){const i=[...n],r=new Array(n.length).fill(0),s=[0,0,0];let o=-1,a=1,l=-1;for(let c=1;c<n.length-1;c++){let m=0;if(n[c]>1e-4){m=Math.atan2(e[c]*o,t[c]*o);let d=0;for(let f=-5;f<=5;f++)Math.abs(m+f*2*Math.PI-s[1])<Math.abs(m+d*2*Math.PI-s[1])&&(d=f);m+=d*2*Math.PI}c>3&&Math.abs(m-s[0])>Math.PI/2&&l<c-1&&n[c]<.5&&(a*=-1,l=c),i[c]*=a,o*=-1,s[0]=s[1],s[1]=m,s[2]=m}const u=[-.086,.343,.486,.343,-.086];for(let c=0;c<n.length-3;c++){let m=0;for(let d=-2;d<=2;d++)m+=i[Math.abs(c+d)]*u[d+2];r[c]=m}for(let c=0;c<n.length-3;c++)i[c]=r[c];const h=7;for(let c=0;c<3;c++){r.fill(0);for(let d=0;d<n.length-h;d++)if(d<h)r[d]=i[d];else{const f=Math.min(5,Math.floor((d-5)/3)),p=ho[f];let g=0;for(let y=-h;y<=h;y++)g+=i[d+y]*p[y+h];r[d]=g}for(let d=n.length-h-2;d<n.length;d++)r[d]=i[d];const m=Math.abs(r[0])>1e-9?r[0]:1;for(let d=0;d<n.length;d++)i[d]=r[d]/m}for(let c=0;c<n.length;c++)i[c]=Math.abs(i[c]);return i}function Go(n,t){const e=[[0,0,0],[0,0,0],[0,0,0]],i=[0,0,0];for(let o=0;o<n.length;o++){const a=n[o],l=-t+o,u=[1,a,a*a];for(let h=0;h<3;h++){i[h]+=u[h]*l;for(let c=0;c<3;c++)e[h][c]+=u[h]*u[c]}}const r=zs(e);if(!r)return null;const s=Ys(r,i);return[s[0],s[1],s[2]]}function Vo(n,t){let e=0,i=1,r=0,s=!1,o=0;const a=Math.min(n.length,se/16*2);for(let l=0;l<a&&!s;l++){const u=n[l];if(i>.5&&u<=.5){const h=-(u-i)*se;Math.abs(h)>1e-9&&(r=-((.5-i-h*e)/h),o=l,s=!0)}i=u,e=l/se}if(!s)return null;if(o>=5&&o<a-10){const l=Math.min(Math.max(2,o-9),9),u=Go(n.slice(o-l,o+l+1),l);if(u){const c=(u[0]+.5*u[1]+.25*u[2]+o)/se;if(o>9)r=c;else{const d=(o-5)/se/8;r=(1-d)*r+d*c}}}return r*Et*t}function Xo(n,t){if(n.length===0)return null;const e=se,i=se/16*4,r=new Fn(e),s=1,o=mo(),a=new Float32Array(501);for(let M=0;M<=500;M++)a[M]=M/500*s*2;const l=new Array(i).fill(0).map((M,v)=>v/e*s*Et),u=new Float32Array(i).fill(0),h=new Float32Array(i).fill(0);let c=0,m=[],d=[],f=[],p=[],g=[],y=[],x=[],b=null,_=0;for(const M of n){const v=M.mtfmapperOrderedDists&&M.mtfmapperOrderedVals&&M.mtfmapperOrderedDists.length===M.mtfmapperOrderedVals.length?bs(M.mtfmapperOrderedDists,M.mtfmapperOrderedVals,M.mtfmapperEffectiveMaxDot??Lt):null,I=(v==null?void 0:v.lsfFull)??M.lsfFull,R=(v==null?void 0:v.esf)??M.esf;if(I.length<e)continue;const U=new Float32Array(e);for(let L=0;L<e;L++)U[L]=I[L]??0;r.transform(U);const Y=Math.max(1e-9,Math.abs(r._real[0])),B=new Array(i).fill(0);for(let L=1;L<i;L++)B[L]=Math.atan2(r._imag[L],r._real[L]);for(let L=0;L<i;L++)u[L]+=r._real[L]/Y,h[L]+=r._imag[L]/Y;if(c++,_+=M.shortSidePx*.5,m.length===0){m=[...I],d=[...R];const L=new Array(i).fill(0);L[0]=1;for(let V=1;V<i;V++)L[V]=Math.hypot(r._real[V]/Y,r._imag[V]/Y);const O=l.map(V=>V),D=(Number.isFinite(M.zeroIndex)?M.zeroIndex:0)*(M.binSize??ui),X=_s(B,O,L,Number.POSITIVE_INFINITY,D),et=Ms(X,l,a);p=et.ptfRaw,g=et.ptfUnwrapped,y=et.ptfLinear,x=et.ptfResidual,f=et.ptfResidual,b=X.fit}}if(c===0)return null;const P=new Float32Array(i),S=new Float32Array(i),w=new Array(i).fill(0);w[0]=1;for(let M=0;M<i;M++)P[M]=u[M]/c,S[M]=h[M]/c,M>0&&(w[M]=Math.hypot(P[M],S[M]));const C=Oo(w,P,S),A=new Array(i).fill(0);for(let M=0;M<i;M++)A[M]=C[M]/o[M];const k=Array.from(a,M=>Ze(M,l,A)),F=Vo(A,s);return{esf:d,lsf:[],lsfCropped:m,mtf:k,ptf:f,ptfRaw:p,ptfUnwrapped:g,ptfLinear:y,ptfResidual:x,ptfPhaseFit:b,freqs:Array.from(a),mtf50:F,calcRadius:_/c}}function zo(n,t,e,i=!1,r=0,s=!1){if(n.length===0)return null;if(n.every(w=>w.mtfmapperLike))return Xo(n);const o=4096,a=new Fn(o),l=1,u=new Float32Array(501);for(let w=0;w<=500;w++)u[w]=w/500*l*2;const h=new Float32Array(501).fill(0);let c=0,m=[],d=[],f=0,p=[],g=[],y=[],x=[],b=[],_=null;for(const w of n){let C=[...w.lsfFull];const A=w.binSize,k=Number.isFinite(w.zeroIndex)?w.zeroIndex:C.length/2,F=Math.max(1,Math.round((w.shortSidePx??0)*.5/Math.max(A,1e-6)));let{peakPos:M,peakIdx:v,peakVal:I}=ni(C,k,F);const R=I*.2;let U=0,Y=C.length-1;for(let T=v;T>=0;T--)if(C[T]<R){U=T;break}for(let T=v;T<C.length;T++)if(C[T]<R){Y=T;break}const B=Y-U;let L=!1;if(i&&B>0){const T=B*4,j=[],G=[];if(r>0){const z=Math.max(0,v-T-r),Z=Math.max(0,v-T);for(let ct=z;ct<Z;ct++)j.push(ct),G.push(C[ct]);const tt=Math.min(C.length,v+T),nt=Math.min(C.length,v+T+r);for(let ct=tt;ct<nt;ct++)j.push(ct),G.push(C[ct])}else{for(let z=0;z<Math.max(0,v-T);z++)j.push(z),G.push(C[z]);for(let z=Math.min(C.length,v+T);z<C.length;z++)j.push(z),G.push(C[z])}if(j.length>2){const{slope:z,intercept:Z}=nr(j,G);for(let tt=0;tt<C.length;tt++)C[tt]=C[tt]-(z*tt+Z);({peakPos:M,peakIdx:v,peakVal:I}=ni(C,k,F)),L=!0}}let O=0,D=0;if(t>0)D=t,O=Math.round(t/A);else{const T=I*.2;let j=0,G=C.length-1;for(let nt=v;nt>=0;nt--)if(C[nt]<T){j=nt;break}for(let nt=v;nt<C.length;nt++)if(C[nt]<T){G=nt;break}const Z=(G-j)*A;let tt=Math.max(2,Z*8);D=tt,O=Math.round(tt/A)}f+=D;const X=Math.max(0,Math.floor(k-O)),et=Math.min(C.length,Math.ceil(k+O)),V=C.slice(X,et);if(V.length<8)continue;const $=new Float32Array(o).fill(0),H=new Array(V.length).fill(0);for(let T=0;T<V.length;T++){let j=1;s&&(j=.5*(1-Math.cos(2*Math.PI*T/(V.length-1)))),H[T]=V[T]*j}const it=Math.max(0,Math.min(V.length-1,M-X));for(let T=0;T<o;T++)$[T]=Co(H,o,T+it);a.transform($);const at=[],E=[],J=[];for(let T=0;T<=o/2;T++){const j=a._real[T],G=a._imag[T],z=Math.sqrt(j*j+G*G);at.push(z),E.push(T/(o*A)*l),J.push(Math.atan2(G,j))}const Q=at[0];if(Q>0){for(let T=0;T<=500;T++){const j=u[T],z=po(j,A);h[T]+=Ze(j,E,at)/Q/z}if(c++,m.length===0){m=Po(V,it,(V.length-1)/2),d=L?Yo(C):w.esf;const T=at.map(Z=>Z/Q),j=E.map(Z=>Z),G=_s(J,j,T,Number.POSITIVE_INFINITY,0),z=Ms(G,E,u);g=z.ptfRaw,y=z.ptfUnwrapped,x=z.ptfLinear,b=z.ptfResidual,p=z.ptfResidual,_=G.fit}}}if(c===0)return null;const P=Array.from(h).map(w=>w/c);let S=null;for(let w=0;w<P.length-1;w++)if(P[w]>=.5&&P[w+1]<.5){S=u[w]+(.5-P[w])*(u[w+1]-u[w])/(P[w+1]-P[w]);break}return{esf:d,lsf:[],lsfCropped:m,mtf:P,ptf:p,ptfRaw:g,ptfUnwrapped:y,ptfLinear:x,ptfResidual:b,ptfPhaseFit:_,freqs:Array.from(u),mtf50:S,calcRadius:f/c}}function Yo(n){const t=new Array(n.length).fill(0);let e=0;for(let i=0;i<n.length;i++)e+=n[i],t[i]=e;return t}function Ze(n,t,e){if(n<=t[0])return e[0];if(n>=t[t.length-1])return e[e.length-1];let i=0;for(;n>t[i+1];)i++;const r=(n-t[i])/(t[i+1]-t[i]);return e[i]+r*(e[i+1]-e[i])}function sr(n){return{...Hi,...n,gradientPercentiles:n!=null&&n.gradientPercentiles&&n.gradientPercentiles.length>0?n.gradientPercentiles:Hi.gradientPercentiles}}function Wo(n){return!n||n.length===0?void 0:[Number.isFinite(n[0])?n[0]:0,Number.isFinite(n[1])?n[1]:Number.isFinite(n[0])?n[0]:0,Number.isFinite(n[2])?n[2]:Number.isFinite(n[0])?n[0]:0,Number.isFinite(n[3])?n[3]:Number.isFinite(n[0])?n[0]:0]}function Ho(n,t){const e=n.width,i=n.height,r=n.data,s=n.bayerPattern||"RGGB",o=Wo(n.blackLevels),a=new Float32Array(e*i),l=(_,P)=>_<0||P<0||_>=e||P>=i?null:Math.max(0,r[P*e+_]-xn(o,_,P));let u=1/0,h=-1/0;for(let _=0;_<i;_++){const P=_*e;for(let S=0;S<e;S++){const w=P+S;let C=0;if(vt(S,_,s,t))C=l(S,_)??0;else{const A=[],k=l(S-1,_),F=l(S+1,_),M=l(S,_-1),v=l(S,_+1);if(k!==null&&vt(S-1,_,s,t)&&A.push(k),F!==null&&vt(S+1,_,s,t)&&A.push(F),M!==null&&vt(S,_-1,s,t)&&A.push(M),v!==null&&vt(S,_+1,s,t)&&A.push(v),A.length>0)C=ei(A);else{const I=[],R=l(S-1,_-1),U=l(S+1,_-1),Y=l(S-1,_+1),B=l(S+1,_+1);R!==null&&vt(S-1,_-1,s,t)&&I.push(R),U!==null&&vt(S+1,_-1,s,t)&&I.push(U),Y!==null&&vt(S-1,_+1,s,t)&&I.push(Y),B!==null&&vt(S+1,_+1,s,t)&&I.push(B),C=ei(I)}}a[w]=C,C<u&&(u=C),C>h&&(h=C)}}if(!Number.isFinite(u)||!Number.isFinite(h)||h<=u+1e-9)return new Uint8Array(e*i);const c=1024,m=new Uint32Array(c),d=h-u;for(let _=0;_<a.length;_++){const P=Math.max(0,Math.min(1,(a[_]-u)/d)),S=Math.min(c-1,Math.max(0,Math.floor(P*(c-1))));m[S]++}const f=a.length,p=_=>{const P=f*_;let S=0;for(let w=0;w<c;w++)if(S+=m[w],S>=P)return u+w/Math.max(1,c-1)*d;return h},g=p(.01),y=p(.99),x=Math.max(1e-9,y-g),b=new Uint8Array(e*i);for(let _=0;_<a.length;_++){const P=Math.max(0,Math.min(1,(a[_]-g)/x));b[_]=Math.round(P*255)}return b}function jo(n,t,e){const i=new Float32Array(n.length),r=new Float32Array(n.length),s=new Float32Array(n.length);for(let o=1;o<e-1;o++)for(let a=1;a<t-1;a++){const l=o*t+a,u=n[(o-1)*t+(a-1)],h=n[(o-1)*t+a],c=n[(o-1)*t+(a+1)],m=n[o*t+(a-1)],d=n[o*t+(a+1)],f=n[(o+1)*t+(a-1)],p=n[(o+1)*t+a],g=n[(o+1)*t+(a+1)],y=-u-2*m-f+(c+2*d+g),x=-u-2*h-c+(f+2*p+g);i[l]=y,r[l]=x,s[l]=Math.hypot(y,x)}return{gx:i,gy:r,magnitude:s}}function Qo(n,t){let e=0,i=0;for(let l=0;l<n.length;l++){const u=n[l];!Number.isFinite(u)||u<=1e-6||(e=Math.max(e,u),i++)}if(i===0||e<=1e-6)return[];const r=1024,s=new Uint32Array(r);for(let l=0;l<n.length;l++){const u=n[l];if(!Number.isFinite(u)||u<=1e-6)continue;const h=Math.max(0,Math.min(1,u/e)),c=Math.min(r-1,Math.floor(h*(r-1)));s[c]++}const o=t&&t.length>0?t:Hi.gradientPercentiles,a=[];for(const l of o){const u=i*l;let h=0;for(let c=0;c<r;c++)if(h+=s[c],h>=u){a.push(c/Math.max(1,r-1)*e);break}}return Array.from(new Set(a.filter(l=>l>0))).sort((l,u)=>u-l)}function qo(n,t){const e=new Uint8Array(n.length);for(let i=0;i<n.length;i++)e[i]=n[i]>=t?1:0;return e}const Ko=256*256;function $o(n,t,e){if(n.length>=Ko){const s=lo.compute(n,t,e);if(s)return{gray:s.blurredGray,gradient:{gx:s.gx,gy:s.gy,magnitude:s.magnitude},backend:"webgl"}}const i=ol(n,t,e),r=jo(i,t,e);return{gray:i,gradient:r,backend:"cpu"}}function Jo(n,t,e,i){let r=n;for(let s=0;s<i;s++){const o=new Uint8Array(n.length);for(let a=0;a<e;a++)for(let l=0;l<t;l++){let u=0;for(let h=-1;h<=1&&!u;h++){const c=a+h;if(!(c<0||c>=e))for(let m=-1;m<=1;m++){const d=l+m;if(!(d<0||d>=t)&&r[c*t+d]){u=1;break}}}o[a*t+l]=u}r=o}return r}function Zo(n,t,e){const i=new Int32Array(n.length),r=[];let s=1;for(let o=0;o<n.length;o++){if(!n[o]||i[o]!==0)continue;const a=[o];i[o]=s;let l=0,u=t,h=e,c=0,m=0,d=0,f=!1;for(;l<a.length;){const p=a[l++],g=p%t,y=Math.floor(p/t);d++,u=Math.min(u,g),h=Math.min(h,y),c=Math.max(c,g),m=Math.max(m,y),(g===0||y===0||g===t-1||y===e-1)&&(f=!0);for(let x=-1;x<=1;x++)for(let b=-1;b<=1;b++){if(b===0&&x===0)continue;const _=g+b,P=y+x;if(_<0||P<0||_>=t||P>=e)continue;const S=P*t+_;!n[S]||i[S]!==0||(i[S]=s,a.push(S))}}r.push({label:s,x:u,y:h,w:c-u+1,h:m-h+1,area:d,touchesBorder:f}),s++}return{labels:i,components:r}}function As(n,t){const e=Math.hypot(n,t);if(!Number.isFinite(e)||e<=1e-9)return null;let i=n/e,r=t/e;return(i<0||Math.abs(i)<=1e-9&&r<0)&&(i=-i,r=-r),{x:i,y:r}}function Ae(n,t){if(n.length===0)return 0;const e=[...n].sort((o,a)=>o.value-a.value),i=e.reduce((o,a)=>o+Math.max(0,a.weight),0);if(i<=0)return e[Math.floor((e.length-1)*t)].value;const r=Math.max(0,Math.min(1,t))*i;let s=0;for(const o of e)if(s+=Math.max(0,o.weight),s>=r)return o.value;return e[e.length-1].value}function Xr(n){const t=n.filter(i=>Number.isFinite(i)).sort((i,r)=>i-r);if(t.length===0)return 0;const e=i=>{if(i.length===1)return i[0];if(i.length===2)return(i[0]+i[1])*.5;const r=Math.ceil(i.length*.5);let s=0,o=1/0;for(let a=0;a+r-1<i.length;a++){const l=i[a+r-1]-i[a];l<o&&(o=l,s=a)}return e(i.slice(s,s+r))};return e(t)}function tl(n,t,e,i,r,s,o){const a=[];for(let l=r.y;l<r.y+r.h;l++)for(let u=r.x;u<r.x+r.w;u++){const h=l*s+u;if(n[h]!==t||!e[h])continue;const c=i.magnitude[h];!Number.isFinite(c)||c<=1e-6||a.push({x:u,y:l,weight:c,gx:i.gx[h],gy:i.gy[h]})}return a}function el(n){let t=0,e=0,i=0,r=0,s=0;for(const l of n){t+=l.weight,e+=l.x*l.weight,i+=l.y*l.weight;const u=Math.hypot(l.gx,l.gy);if(!Number.isFinite(u)||u<=1e-6)continue;const h=-l.gy/u,c=l.gx/u;r+=l.weight*(h*h-c*c),s+=l.weight*(2*h*c)}if(t<=0)return null;e/=t,i/=t;const o=.5*Math.atan2(s,r),a=As(Math.cos(o),Math.sin(o));return a?{centerX:e,centerY:i,dirX:a.x,dirY:a.y,orthoX:-a.y,orthoY:a.x}:null}function Wn(n,t){let e=0,i=0;const r=-t.dirY,s=t.dirX;for(const o of n){const a=(o.x-t.pointX)*r+(o.y-t.pointY)*s;i+=o.weight*a*a,e+=o.weight}return e<=0?1/0:Math.sqrt(i/e)}function zr(n,t,e,i,r){const s=Math.max(0,Math.min(t-1,i)),o=Math.max(0,Math.min(e-1,r)),a=Math.floor(s),l=Math.floor(o),u=Math.min(t-1,a+1),h=Math.min(e-1,l+1),c=s-a,m=o-l,d=n[l*t+a],f=n[l*t+u],p=n[h*t+a],g=n[h*t+u],y=d+(f-d)*c,x=p+(g-p)*c;return y+(x-y)*m}function nl(n,t,e,i,r,s,o,a){const l=zr(n,t,e,i-s*a,r-o*a);return zr(n,t,e,i+s*a,r+o*a)-l}function Hn(n,t,e,i,r){const s=Math.max(1e-6,e-t);if(n.length===0||!Number.isFinite(s))return{points:[],coverageRatio:0,centerCoverageRatio:0};const o=Math.max(1.5,Math.min(4,s/18)),a=Math.max(1,Math.ceil(s/o)),l=new Map;for(const f of n){const p=i(f);if(!Number.isFinite(p)||p<t||p>e)continue;const g=Math.max(0,Math.min(a-1,Math.floor((p-t)/o))),y=f.weight/(1+Math.abs(r(f))),x=l.get(g);(!x||y>x.score)&&l.set(g,{point:f,score:y})}const u=Array.from(l.values()).sort((f,p)=>i(f.point)-i(p.point)).map(f=>f.point),h=Math.max(0,Math.floor(a*.3)),c=Math.max(h+1,Math.ceil(a*.7));let m=0;for(let f=h;f<c;f++)l.has(f)&&m++;const d=Math.max(1,c-h);return{points:u,coverageRatio:u.length/a,centerCoverageRatio:m/d}}function jn(n,t){const e=n.dirX*t.dirY-n.dirY*t.dirX;if(!Number.isFinite(e)||Math.abs(e)<=1e-6)return null;const i=t.pointX-n.pointX,r=t.pointY-n.pointY,s=(i*t.dirY-r*t.dirX)/e;return{x:n.pointX+n.dirX*s,y:n.pointY+n.dirY*s}}function il(n){if(n.length<3)return 0;let t=0;for(let e=0;e<n.length;e++){const i=n[e],r=n[(e+1)%n.length];t+=i.x*r.y-r.x*i.y}return t*.5}function rl(n,t,e,i,r,s,o,a,l){const u=tl(i,r.label,s,o,r,t),h=u.map(N=>({x:N.x,y:N.y}));if(u.length<l.minEdgePoints)return{candidate:null,failureStage:"min_edge_points",pointsCount:u.length,strongEdgePoints:h};const c=el(u);if(!c)return{candidate:null,failureStage:"dominant_axes",pointsCount:u.length,strongEdgePoints:h};const m=u.map(N=>{const lt=N.x-c.centerX,ht=N.y-c.centerY;return{...N,u:lt*c.dirX+ht*c.dirY,v:lt*c.orthoX+ht*c.orthoY}}),d={x:c.centerX,y:c.centerY},f=Ae(m.map(N=>({value:N.u,weight:N.weight})),l.extentQuantileLow),p=Ae(m.map(N=>({value:N.u,weight:N.weight})),l.extentQuantileHigh),g=Ae(m.map(N=>({value:N.v,weight:N.weight})),l.extentQuantileLow),y=Ae(m.map(N=>({value:N.v,weight:N.weight})),l.extentQuantileHigh),x=Math.max(1e-6,Math.max(Math.abs(f),Math.abs(p))),b=Math.max(1e-6,Math.max(Math.abs(g),Math.abs(y))),_=72,P=360/_,S=Array.from({length:_},()=>[]),w=N=>{let lt=N%360;return lt<0&&(lt+=360),lt},C=(N,lt)=>{const ht=Math.abs(w(N)-w(lt));return Math.min(ht,360-ht)};m.forEach(N=>{const lt=N.u/x,ht=N.v/b,At=w(Math.atan2(ht,lt)*180/Math.PI),Yt=Math.hypot(lt,ht),Rt=Math.max(0,Math.min(_-1,Math.floor(At/P)));S[Rt].push({point:N,angleDeg:At,normRadius:Yt})});const A=S.map(N=>N.length>0?Xr(N.map(lt=>lt.normRadius)):-1/0),k=(N,lt)=>{let ht=-1,At=-1/0;for(let yt=0;yt<S.length;yt++){if(S[yt].length===0)continue;const pe=(yt+.5)*P;if(C(pe,N)>45||lt.some(Gs=>C(pe,Gs)<45))continue;const Fe=A[yt];Fe>At&&(At=Fe,ht=yt)}let Yt=ht>=0?(ht+.5)*P:N,Rt=ht>=0?S[ht]:m.map(yt=>{const fe=yt.u/x,pe=yt.v/b;return{point:yt,angleDeg:w(Math.atan2(pe,fe)*180/Math.PI),normRadius:Math.hypot(fe,pe)}}).filter(yt=>C(yt.angleDeg,N)<=45&&!lt.some(fe=>C(yt.angleDeg,fe)<45));if(Rt.length===0&&(Rt=m.map(yt=>{const fe=yt.u/x,pe=yt.v/b;return{point:yt,angleDeg:w(Math.atan2(pe,fe)*180/Math.PI),normRadius:Math.hypot(fe,pe)}}).filter(yt=>C(yt.angleDeg,N)<=45),Yt=N),Rt.length===0)return{x:m[0].x,y:m[0].y,u:m[0].u,v:m[0].v,angleDeg:N};const Gn=ht>=0?A[ht]:Xr(Rt.map(yt=>yt.normRadius));let Me=0,wn=0,pr=0,mr=0,gr=0;for(const yt of Rt){const fe=C(yt.angleDeg,N)/45,pe=Math.abs(yt.normRadius-Gn),Fe=Math.max(1e-6,yt.point.weight)/(1+fe*2+pe*6);Me+=Fe,wn+=yt.point.x*Fe,pr+=yt.point.y*Fe,mr+=yt.point.u*Fe,gr+=yt.point.v*Fe}return Me>0?{x:wn/Me,y:pr/Me,u:mr/Me,v:gr/Me,angleDeg:Yt}:{x:Rt[0].point.x,y:Rt[0].point.y,u:Rt[0].point.u,v:Rt[0].point.v,angleDeg:Rt[0].angleDeg}},F=k(225,[]),M=k(315,[F.angleDeg]),v=k(45,[F.angleDeg,M.angleDeg]),I=k(135,[F.angleDeg,M.angleDeg,v.angleDeg]),R=[{x:F.x,y:F.y},{x:M.x,y:M.y},{x:v.x,y:v.y},{x:I.x,y:I.y}],U=p-f,Y=y-g,B=Math.min(U,Y),L=Math.max(U,Y);if(!Number.isFinite(B)||B<l.minSpanPx||L/Math.max(1,B)>l.maxAspectRatio)return{candidate:null,failureStage:"span_aspect",pointsCount:u.length,minSpan:B,maxSpan:L,axisCentroid:d,axisExtremePoints:R,strongEdgePoints:h};const O=Math.max(l.bandMinPx,Math.min(l.bandMaxPx,B*l.bandScale)),D=Math.max(1,Math.min(3,O*.55)),X=Math.max(a,0),et=void 0,V=void 0,$=N=>N.map(lt=>({x:lt.x,y:lt.y,weight:lt.weight})),H=N=>N.map(lt=>({x:lt.x,y:lt.y})),it=(N,lt,ht)=>N.filter(At=>{if(!Number.isFinite(At.weight)||At.weight<X)return!1;const Yt=nl(n,t,e,At.x,At.y,lt,ht,D);return Number.isFinite(Yt)&&Yt>=l.minPointContrast}),at=f,E=p,J=g,Q=y,T=l.minCoverageRatio,j=l.minCenterCoverageRatio,G=[],z=[],Z=[],tt=[],nt=[],ct=(N,lt,ht,At,Yt,Rt)=>(ht-N)*(Rt-lt)-(At-lt)*(Yt-N),dt=N=>N>1e-6?1:N<-1e-6?-1:0,rt=[{u:(F.u+M.u)*.5,v:(F.v+M.v)*.5},{u:(M.u+v.u)*.5,v:(M.v+v.v)*.5},{u:(v.u+I.u)*.5,v:(v.v+I.v)*.5},{u:(I.u+F.u)*.5,v:(I.v+F.v)*.5}],Bt=(N,lt)=>{const ht=dt(ct(F.u,F.v,v.u,v.v,N,lt)),At=dt(ct(M.u,M.v,I.u,I.v,N,lt));return`${ht},${At}`},Vt=new Map;rt.forEach((N,lt)=>{Vt.set(Bt(N.u,N.v),lt)});for(const N of m){if(!Number.isFinite(N.u)||!Number.isFinite(N.v)){nt.push(N);continue}let ht=Vt.get(Bt(N.u,N.v))??-1;if(ht<0){let At=1/0;for(let Yt=0;Yt<rt.length;Yt++){const Rt=rt[Yt],Gn=(N.u-Rt.u)/x,Me=(N.v-Rt.v)/b,wn=Gn*Gn+Me*Me;wn<At&&(At=wn,ht=Yt)}}ht===0?G.push(N):ht===1?z.push(N):ht===2?Z.push(N):ht===3?tt.push(N):nt.push(N)}const W=[...G,...Z],ut=[...z,...tt],ot={dir:W.length,ortho:ut.length,unassigned:m.length-W.length-ut.length},pt=G.length>=l.minSidePoints?Ae(G.map(N=>({value:N.v,weight:N.weight})),.5):g,mt=Z.length>=l.minSidePoints?Ae(Z.map(N=>({value:N.v,weight:N.weight})),.5):y,te=tt.length>=l.minSidePoints?Ae(tt.map(N=>({value:N.u,weight:N.weight})),.5):f,Dt=z.length>=l.minSidePoints?Ae(z.map(N=>({value:N.u,weight:N.weight})),.5):p,_t=[{x:(F.x+M.x)*.5,y:(F.y+M.y)*.5},{x:(M.x+v.x)*.5,y:(M.y+v.y)*.5},{x:(v.x+I.x)*.5,y:(v.y+I.y)*.5},{x:(I.x+F.x)*.5,y:(I.y+F.y)*.5}],le=G.filter(N=>Math.abs(N.v-pt)<=O&&N.u>=at&&N.u<=E),Kt=Z.filter(N=>Math.abs(N.v-mt)<=O&&N.u>=at&&N.u<=E),ke=tt.filter(N=>Math.abs(N.u-te)<=O&&N.v>=J&&N.v<=Q),K=z.filter(N=>Math.abs(N.u-Dt)<=O&&N.v>=J&&N.v<=Q),Ft=[le.length,K.length,Kt.length,ke.length],Tt=[H(le),H(K),H(Kt),H(ke)],xt=it(le,-c.orthoX,-c.orthoY),wt=it(Kt,c.orthoX,c.orthoY),ee=it(ke,-c.dirX,-c.dirY),$t=it(K,c.dirX,c.dirY),_n=[xt.length,$t.length,wt.length,ee.length],Be=[H(xt),H($t),H(wt),H(ee)],De=Hn(xt,at,E,N=>N.u,N=>N.v-pt),Oe=Hn($t,J,Q,N=>N.v,N=>N.u-Dt),Ge=Hn(wt,at,E,N=>N.u,N=>N.v-mt),Ve=Hn(ee,J,Q,N=>N.v,N=>N.u-te),Tn=(N,lt)=>N.slice().sort((ht,At)=>lt(ht)-lt(At)),In=Tn(xt,N=>N.u),Nn=Tn($t,N=>N.v),Rn=Tn(wt,N=>N.u),Ln=Tn(ee,N=>N.v),di=[In.length,Nn.length,Rn.length,Ln.length],Ht=[De.coverageRatio,Oe.coverageRatio,Ge.coverageRatio,Ve.coverageRatio];De.centerCoverageRatio,Oe.centerCoverageRatio,Ge.centerCoverageRatio,Ve.centerCoverageRatio;const Xt=[H(In),H(Nn),H(Rn),H(Ln)],Jt={axisPointCounts:ot,sideBandPointCounts:Ft,sideContrastPointCounts:_n,gradientThreshold:a,pointAxisMinDot:et,pointAxisMargin:V,bandWidth:O,minPointContrast:l.minPointContrast,minCoverageRatio:T,minCenterCoverageRatio:j,axisCentroid:d,axisExtremePoints:R,axisSideCenters:_t,strongEdgePoints:h,axisDirPoints:H(W),axisOrthoPoints:H(ut),axisUnassignedPoints:H(nt),sideBandPoints:Tt,sideContrastPoints:Be};if(In.length<l.minSidePoints||Rn.length<l.minSidePoints||Ln.length<l.minSidePoints||Nn.length<l.minSidePoints)return{candidate:null,failureStage:"min_side_points",pointsCount:u.length,minSpan:B,maxSpan:L,sidePointCounts:di,sideCoverageRatios:Ht,...Jt,sideFitPoints:Xt};if(De.coverageRatio<T||Oe.coverageRatio<T||Ge.coverageRatio<T||Ve.coverageRatio<T||De.centerCoverageRatio<j||Oe.centerCoverageRatio<j||Ge.centerCoverageRatio<j||Ve.centerCoverageRatio<j)return{candidate:null,failureStage:"side_coverage",pointsCount:u.length,minSpan:B,maxSpan:L,sidePointCounts:di,sideCoverageRatios:Ht,...Jt,sideFitPoints:Xt};const fi=In,pi=Rn,mi=Ln,gi=Nn,ne=di,Xe=ve($(fi)),ze=ve($(pi)),Ye=ve($(mi)),We=ve($(gi));if(!Xe||!ze||!Ye||!We)return{candidate:null,failureStage:"fit_lines",pointsCount:u.length,minSpan:B,maxSpan:L,sidePointCounts:ne,sideCoverageRatios:Ht,...Jt,sideFitPoints:Xt};const _e=Lo($(fi),Xe,$(gi),We,$(pi),ze,$(mi),Ye),En=l.minAxisDot,Un=(N,lt,ht)=>Math.abs(N.dirX*lt+N.dirY*ht),zt=[Un(Xe,c.dirX,c.dirY),Un(We,c.orthoX,c.orthoY),Un(ze,c.dirX,c.dirY),Un(Ye,c.orthoX,c.orthoY)];if(zt[0]<En||zt[1]<En||zt[2]<En||zt[3]<En)return{candidate:null,failureStage:"axis_alignment",pointsCount:u.length,minSpan:B,maxSpan:L,sidePointCounts:ne,sideCoverageRatios:Ht,axisDots:zt,...Jt,sideFitPoints:Xt,sideFitLines:_e};const ie=Math.max(l.residualLimitFloor,O*l.residualLimitScale),gt=[Wn($(fi),Xe),Wn($(pi),ze),Wn($(mi),Ye),Wn($(gi),We)],yi=[gt[0],gt[3],gt[1],gt[2]],xi=Math.max(...gt),ce=jn(Xe,Ye),ue=jn(Xe,We),he=jn(ze,We),de=jn(ze,Ye);if(!ce||!ue||!he||!de)return{candidate:null,failureStage:"corners",pointsCount:u.length,minSpan:B,maxSpan:L,sidePointCounts:ne,sideCoverageRatios:Ht,axisDots:zt,sideResiduals:[gt[0],gt[3],gt[1],gt[2]],residualLimit:ie,...Jt,sideFitPoints:Xt,sideFitLines:_e};const Mn=[ce,ue,he,de],re=Math.abs(il(Mn));if(!Number.isFinite(re)||re<l.minQuadArea)return{candidate:null,failureStage:"quad_area",pointsCount:u.length,minSpan:B,maxSpan:L,sidePointCounts:ne,sideCoverageRatios:Ht,axisDots:zt,sideResiduals:[gt[0],gt[3],gt[1],gt[2]],residualLimit:ie,...Jt,sideFitPoints:Xt,sideFitLines:_e,quadArea:re};const bi=Math.hypot(ue.x-ce.x,ue.y-ce.y),_i=Math.hypot(he.x-ue.x,he.y-ue.y),Mi=Math.hypot(he.x-de.x,he.y-de.y),wi=Math.hypot(de.x-ce.x,de.y-ce.y),He=[bi,_i,Mi,wi],Si=Math.min(bi,_i,Mi,wi),Us=Math.max(bi,_i,Mi,wi);if(!Number.isFinite(Si)||Si<l.minSideLength||Us/Math.max(1,Si)>l.maxAspectRatio)return{candidate:null,failureStage:"side_length",pointsCount:u.length,minSpan:B,maxSpan:L,sidePointCounts:ne,sideCoverageRatios:Ht,axisDots:zt,sideResiduals:[gt[0],gt[3],gt[1],gt[2]],residualLimit:ie,...Jt,sideFitPoints:Xt,sideFitLines:_e,quadArea:re,sideLengths:He};const je=As(ue.x-ce.x+(he.x-de.x),ue.y-ce.y+(he.y-de.y));if(!je)return{candidate:null,failureStage:"corners",pointsCount:u.length,minSpan:B,maxSpan:L,sidePointCounts:ne,sideCoverageRatios:Ht,axisDots:zt,sideResiduals:[gt[0],gt[3],gt[1],gt[2]],residualLimit:ie,...Jt,sideFitPoints:Xt,sideFitLines:_e,quadArea:re,sideLengths:He};const hr={x:-je.y,y:je.x},vi=(ce.x+ue.x+he.x+de.x)*.25,Ci=(ce.y+ue.y+he.y+de.y)*.25,Bn=Mn.map(N=>{const lt=N.x-vi,ht=N.y-Ci;return{u:lt*je.x+ht*je.y,v:lt*hr.x+ht*hr.y}}),Dn=(Math.max(...Bn.map(N=>N.u))-Math.min(...Bn.map(N=>N.u)))*.5,On=(Math.max(...Bn.map(N=>N.v))-Math.min(...Bn.map(N=>N.v)))*.5;if(!Number.isFinite(Dn)||!Number.isFinite(On)||Math.min(Dn,On)<6)return{candidate:null,failureStage:"box_size",pointsCount:u.length,minSpan:B,maxSpan:L,sidePointCounts:ne,sideCoverageRatios:Ht,axisDots:zt,sideResiduals:[gt[0],gt[3],gt[1],gt[2]],residualLimit:ie,sideFitPoints:Xt,quadArea:re,sideLengths:He};const q=hl(n,t,e,Mn,Xe,We,ze,Ye,vi,Ci,Dn,On,l.innerPurityStdScale,l.outerMeanSpreadLimit);if(!Number.isFinite(xi)||xi>ie)return{candidate:null,failureStage:"residual",pointsCount:u.length,minSpan:B,maxSpan:L,sidePointCounts:ne,sideCoverageRatios:Ht,axisDots:zt,sideResiduals:yi,residualLimit:ie,...Jt,sideFitPoints:Xt,sideFitLines:_e,quadArea:re,sideLengths:He,outerContrast:q.contrast,outerUniformityOk:q.ok,outerMeanSpread:q.meanSpread,outerMeanSpreadLimit:q.meanSpreadLimit,outerAvgStd:q.avgStd,outerAvgStdLimit:q.avgStdLimit,outerSideMeans:q.outerSideMeans,outerSideStds:q.outerSideStds,outerSideStdLimit:q.outerSideStdLimit,outerSideQuads:q.outerSideQuads,innerSideUniformityOk:q.innerSideOk,innerSideStds:q.innerSideStds,innerSideStdLimit:q.innerSideStdLimit,innerSideQuads:q.innerSideQuads};const dr=l.filterBlockPurity&&(!q.ok||!q.innerSideOk);if(dr||q.contrast<l.minOuterContrast)return{candidate:null,failureStage:dr?q.ok?"inner_roi_uniformity":"outer_uniformity":"outer_contrast",pointsCount:u.length,minSpan:B,maxSpan:L,sidePointCounts:ne,sideCoverageRatios:Ht,axisDots:zt,sideResiduals:yi,residualLimit:ie,...Jt,sideFitPoints:Xt,sideFitLines:_e,quadArea:re,sideLengths:He,outerContrast:q.contrast,outerUniformityOk:q.ok,outerMeanSpread:q.meanSpread,outerMeanSpreadLimit:q.meanSpreadLimit,outerAvgStd:q.avgStd,outerAvgStdLimit:q.avgStdLimit,outerSideMeans:q.outerSideMeans,outerSideStds:q.outerSideStds,outerSideStdLimit:q.outerSideStdLimit,outerSideQuads:q.outerSideQuads,innerSideUniformityOk:q.innerSideOk,innerSideStds:q.innerSideStds,innerSideStdLimit:q.innerSideStdLimit,innerSideQuads:q.innerSideQuads};const fr=Mt(Gt(Mn,1),t,e);if(!fr)return{candidate:null,failureStage:"bbox",pointsCount:u.length,minSpan:B,maxSpan:L,sidePointCounts:ne,sideCoverageRatios:Ht,axisDots:zt,sideResiduals:[gt[0],gt[3],gt[1],gt[2]],residualLimit:ie,...Jt,sideFitPoints:Xt,sideFitLines:_e,quadArea:re,sideLengths:He,outerContrast:q.contrast,outerUniformityOk:q.ok,outerMeanSpread:q.meanSpread,outerMeanSpreadLimit:q.meanSpreadLimit,outerAvgStd:q.avgStd,outerAvgStdLimit:q.avgStdLimit,outerSideMeans:q.outerSideMeans,outerSideStds:q.outerSideStds,outerSideStdLimit:q.outerSideStdLimit,outerSideQuads:q.outerSideQuads,innerSideUniformityOk:q.innerSideOk,innerSideStds:q.innerSideStds,innerSideStdLimit:q.innerSideStdLimit,innerSideQuads:q.innerSideQuads};const Bs=1/(1+xi/Math.max(1,ie)),Ds=l.filterBlockPurity?q.score:1,Os=q.contrast*Ds*Bs*Math.sqrt(re);return{candidate:{centerX:vi,centerY:Ci,dirX:je.x,dirY:je.y,halfWidth:Dn,halfHeight:On,score:Os,bbox:fr,corners:Mn,sideFitPoints:Xt,outerSideMeans:q.outerSideMeans,outerSideQuads:q.outerSideQuads},failureStage:null,pointsCount:u.length,minSpan:B,maxSpan:L,sidePointCounts:ne,sideCoverageRatios:Ht,axisDots:zt,sideResiduals:yi,residualLimit:ie,...Jt,sideFitPoints:Xt,sideFitLines:_e,quadArea:re,sideLengths:He,outerContrast:q.contrast,outerUniformityOk:q.ok,outerMeanSpread:q.meanSpread,outerMeanSpreadLimit:q.meanSpreadLimit,outerAvgStd:q.avgStd,outerAvgStdLimit:q.avgStdLimit,outerSideMeans:q.outerSideMeans,outerSideStds:q.outerSideStds,outerSideStdLimit:q.outerSideStdLimit,outerSideQuads:q.outerSideQuads,innerSideUniformityOk:q.innerSideOk,innerSideStds:q.innerSideStds,innerSideStdLimit:q.innerSideStdLimit,innerSideQuads:q.innerSideQuads}}function sl(n,t,e,i,r,s,o,a,l){return rl(n,t,e,i,r,s,o,a,l).candidate}function al(n,t,e,i,r,s){const o=sr(r),a=Math.max(i*8,i+128),l=ll(n,t,e,o.downsampleMaxSide);s==null||s("Detecting candidates: downsampling...",.02),s==null||s("Detecting candidates: edge stage...",.06);const u=$o(l.gray,l.width,l.height),h=u.gray,c=u.gradient;s==null||s(`Detecting candidates: gradient (${u.backend==="webgl"?"WebGL1":"CPU"})...`,.1);const m=Qo(c.magnitude,o.gradientPercentiles),d=l.width*l.height,f=Math.max(o.minComponentAreaPx,Math.round(d*o.minComponentAreaRatio)),p=Math.max(f+1,Math.round(d*o.maxComponentAreaRatio)),g=[],y=Math.max(1,m.reduce((F,M,v)=>F+(v<=1,2),0));let x=0;for(let F=0;F<m.length;F++){const M=m[F],v=qo(c.magnitude,M),I=F<=1?[3,2]:[2,1];for(const R of I){const U=x/y;s==null||s(`Detecting candidates: threshold ${F+1}/${m.length}, dilate ${R}`,.12+.78*U);const Y=Jo(v,l.width,l.height,R),{labels:B,components:L}=Zo(Y,l.width,l.height);for(const O of L){if(O.touchesBorder||O.area<f||O.area>p)continue;const D=sl(h,l.width,l.height,B,O,v,c,M,o);if(!D)continue;const X=1/l.scale,et=D.corners.map(V=>({x:V.x*X,y:V.y*X}));g.push({centerX:D.centerX*X,centerY:D.centerY*X,dirX:D.dirX,dirY:D.dirY,halfWidth:D.halfWidth*X,halfHeight:D.halfHeight*X,score:D.score,bbox:{x:D.bbox.x*X,y:D.bbox.y*X,w:D.bbox.w*X,h:D.bbox.h*X},corners:et,sideFitPoints:D.sideFitPoints?[D.sideFitPoints[0].map(V=>({x:V.x*X,y:V.y*X})),D.sideFitPoints[1].map(V=>({x:V.x*X,y:V.y*X})),D.sideFitPoints[2].map(V=>({x:V.x*X,y:V.y*X})),D.sideFitPoints[3].map(V=>({x:V.x*X,y:V.y*X}))]:void 0,outerSideMeans:D.outerSideMeans,outerSideQuads:D.outerSideQuads?[D.outerSideQuads[0].map(V=>({x:V.x*X,y:V.y*X})),D.outerSideQuads[1].map(V=>({x:V.x*X,y:V.y*X})),D.outerSideQuads[2].map(V=>({x:V.x*X,y:V.y*X})),D.outerSideQuads[3].map(V=>({x:V.x*X,y:V.y*X}))]:void 0})}g.length>a&&(g.sort((O,D)=>D.score-O.score),g.length=a),x++}}console.log(`[SFR Auto Detect] Candidate pool before dedupe: ${g.length}`),s==null||s(`Detecting candidates: deduplicating (0/${Math.max(1,Math.min(g.length,Math.max(i*4,i+32)))})...`,.94),g.sort((F,M)=>M.score-F.score);const b=Math.max(i*4,i+32),_=g.length>b?g.slice(0,b):g,P=[];if(_.length<=256){console.log(`[SFR Auto Detect] Using simple dedupe for ${_.length} candidates`);for(let F=0;F<_.length;F++){const M=_[F];console.log(`[SFR Auto Detect] Simple dedupe candidate ${F+1}/${_.length}`,M.bbox);const v=_.length<=0?1:F/_.length;if(s==null||s(`Detecting candidates: deduplicating (${F}/${_.length})...`,.94+.05*Math.min(1,v)),!P.some(R=>{const U=Math.hypot(M.centerX-R.centerX,M.centerY-R.centerY),Y=Math.max(Math.hypot(M.bbox.w,M.bbox.h),Math.hypot(R.bbox.w,R.bbox.h));return Yr(M.bbox,R.bbox)>.28||U<Y*.18})&&(P.push(M),P.length>=i))break}return s==null||s("Detecting candidates: deduplicating...",1),P}const S=Math.max(32,Math.round(Math.sqrt(Math.max(1,t*e)/4096))),w=new Map,C=new Set,A=F=>Math.floor(F/S),k=(F,M)=>{if(!Number.isFinite(F.bbox.x)||!Number.isFinite(F.bbox.y)||!Number.isFinite(F.bbox.w)||!Number.isFinite(F.bbox.h)||F.bbox.w<=0||F.bbox.h<=0||F.bbox.w>t*4||F.bbox.h>e*4)return;const v=A(F.bbox.x),I=A(F.bbox.x+F.bbox.w),R=A(F.bbox.y),U=A(F.bbox.y+F.bbox.h);for(let Y=R;Y<=U;Y++)for(let B=v;B<=I;B++){const L=`${B},${Y}`,O=w.get(L);O?O.push(M):w.set(L,[M])}};for(let F=0;F<_.length;F++){const M=_[F];if(F===0||F%200===0){const B=_.length<=0?1:F/_.length;s==null||s(`Detecting candidates: deduplicating (${F}/${_.length})...`,.94+.05*Math.min(1,B))}C.clear();const v=A(M.bbox.x),I=A(M.bbox.x+M.bbox.w),R=A(M.bbox.y),U=A(M.bbox.y+M.bbox.h);let Y=!1;for(let B=R-1;B<=U+1&&!Y;B++)for(let L=v-1;L<=I+1&&!Y;L++){const O=w.get(`${L},${B}`);if(O)for(const D of O){if(C.has(D))continue;C.add(D);const X=P[D];if(!X)continue;const et=Math.hypot(M.centerX-X.centerX,M.centerY-X.centerY),V=Math.max(Math.hypot(M.bbox.w,M.bbox.h),Math.hypot(X.bbox.w,X.bbox.h));if(Yr(M.bbox,X.bbox)>.28||et<V*.18){Y=!0;break}}}if(!Y){const B=P.length;if(P.push(M),k(M,B),P.length>=i)break}}return s==null||s("Detecting candidates: deduplicating...",1),P}function ol(n,t,e){const i=new Uint8Array(n.length);for(let r=0;r<e;r++)for(let s=0;s<t;s++){let o=0,a=0;for(let l=-1;l<=1;l++){const u=r+l;if(!(u<0||u>=e))for(let h=-1;h<=1;h++){const c=s+h;c<0||c>=t||(o+=n[u*t+c],a++)}}i[r*t+s]=Math.round(o/Math.max(1,a))}return i}function ll(n,t,e,i){const r=Math.max(t,e);if(r<=i)return{gray:n,width:t,height:e,scale:1};const s=i/r,o=Math.max(1,Math.round(t*s)),a=Math.max(1,Math.round(e*s)),l=new Uint8Array(o*a);for(let u=0;u<a;u++){const h=Math.min(e-1,Math.floor(u/s));for(let c=0;c<o;c++){const m=Math.min(t-1,Math.floor(c/s));l[u*o+c]=n[h*t+m]}}return{gray:l,width:o,height:a,scale:s}}function Yr(n,t){const e=Math.max(n.x,t.x),i=Math.max(n.y,t.y),r=Math.min(n.x+n.w,t.x+t.w),s=Math.min(n.y+n.h,t.y+t.h),o=Math.max(0,r-e),a=Math.max(0,s-i),l=o*a;if(l<=0)return 0;const u=n.w*n.h+t.w*t.h-l;return u>0?l/u:0}function Wr(n){const t=n.length;if(t===0)return{count:0,mean:0,std:1/0};let e=0;for(const s of n)e+=s;const i=e/t;let r=0;for(const s of n){const o=s-i;r+=o*o}return r/=t,{count:t,mean:i,std:Math.sqrt(Math.max(0,r))}}function cl(n,t,e,i){return{p1:{x:n.x-t*i,y:n.y-e*i},p2:{x:n.x+t*i,y:n.y+e*i}}}function Hr(n,t,e,i,r){return[{x:n.p1.x+t*i,y:n.p1.y+e*i},{x:n.p2.x+t*i,y:n.p2.y+e*i},{x:n.p2.x+t*r,y:n.p2.y+e*r},{x:n.p1.x+t*r,y:n.p1.y+e*r}]}function ul(n,t,e){let i=0;for(let r=0;r<4;r++){const s=e[r],o=e[(r+1)%4],a=(o.x-s.x)*(t-s.y)-(o.y-s.y)*(n-s.x);if(Math.abs(a)<=1e-6)continue;const l=a>0?1:-1;if(i===0)i=l;else if(i!==l)return!1}return!0}function jr(n,t,e,i){const r=Mt(Gt(i,1),t,e);if(!r)return[];const s=[];for(let o=r.y;o<r.y+r.h;o++)for(let a=r.x;a<r.x+r.w;a++)ul(a,o,i)&&s.push(n[o*t+a]);return s}function hl(n,t,e,i,r,s,o,a,l,u,h,c,m,d){const f=h*2,p=c*2,g=Math.hypot(i[1].x-i[0].x,i[1].y-i[0].y),y=Math.hypot(i[2].x-i[1].x,i[2].y-i[1].y),x=Math.hypot(i[2].x-i[3].x,i[2].y-i[3].y),b=Math.hypot(i[3].x-i[0].x,i[3].y-i[0].y),P=Math.max(...[g,y,x,b]),S=Math.max(2,Math.min(f,p)),w=Math.max(4,P*.25),C=Math.max(2,Math.min(12,S*.22)),A=Math.max(1,Math.min(C,Math.max(1,S*.5-1))),k=1,F=Math.max(8,Math.round(Math.min(w,C*3))),M=[[i[0],i[1],i[1],i[0]],[i[1],i[2],i[2],i[1]],[i[2],i[3],i[3],i[2]],[i[3],i[0],i[0],i[3]]],v=[[i[0],i[1],i[1],i[0]],[i[1],i[2],i[2],i[1]],[i[2],i[3],i[3],i[2]],[i[3],i[0],i[0],i[3]]],I=[],R=[],U=[{corners:[i[0],i[1]],seedLine:r,sideLength:g},{corners:[i[1],i[2]],seedLine:s,sideLength:y},{corners:[i[2],i[3]],seedLine:o,sideLength:x},{corners:[i[3],i[0]],seedLine:a,sideLength:b}];for(let G=0;G<U.length;G++){const z=U[G],Z=Math.max(1,z.sideLength*.5-1),tt=Math.max(1,Math.min(Z,w*.5)),nt={x:(z.corners[0].x+z.corners[1].x)*.5,y:(z.corners[0].y+z.corners[1].y)*.5},ct=cl(nt,z.seedLine.dirX,z.seedLine.dirY,tt),dt=ge(n,t,e,ct.p1,ct.p2,tt,Math.max(4,k+Math.max(C,A)+2)),rt=(dt==null?void 0:dt.line)||ct,Bt=rt.p2.x-rt.p1.x,Vt=rt.p2.y-rt.p1.y,W=Math.hypot(Bt,Vt);if(!Number.isFinite(W)||W<=1e-6)return{ok:!1,score:0,meanSpread:1/0,meanSpreadLimit:1/0,avgStd:1/0,avgStdLimit:1/0,contrast:0,outerMean:0,outerSideMeans:[0,0,0,0],outerSideStds:[1/0,1/0,1/0,1/0],outerSideStdLimit:1/0,outerSideQuads:M,innerSideOk:!1,innerSideStds:[1/0,1/0,1/0,1/0],innerSideStdLimit:1/0,innerSideQuads:v};let ut=-Vt/W,ot=Bt/W;const pt={x:(rt.p1.x+rt.p2.x)*.5,y:(rt.p1.y+rt.p2.y)*.5};(pt.x-l)*ut+(pt.y-u)*ot<0&&(ut=-ut,ot=-ot);const mt=Hr(rt,ut,ot,k,k+C),te=Hr(rt,ut,ot,-k,-(k+A));M[G]=mt,v[G]=te,I.push(jr(n,t,e,mt)),R.push(jr(n,t,e,te))}const Y=I.map(Wr);if(Y.some(G=>G.count<F||!Number.isFinite(G.std)))return{ok:!1,score:0,meanSpread:1/0,meanSpreadLimit:1/0,avgStd:1/0,avgStdLimit:1/0,contrast:0,outerMean:0,outerSideMeans:[0,0,0,0],outerSideStds:[1/0,1/0,1/0,1/0],outerSideStdLimit:1/0,outerSideQuads:M,innerSideOk:!1,innerSideStds:[1/0,1/0,1/0,1/0],innerSideStdLimit:1/0,innerSideQuads:v};const B=R.map(Wr);if(B.some(G=>G.count<F||!Number.isFinite(G.std)||!Number.isFinite(G.mean)))return{ok:!1,score:0,meanSpread:1/0,meanSpreadLimit:1/0,avgStd:1/0,avgStdLimit:1/0,contrast:0,outerMean:0,outerSideMeans:[0,0,0,0],outerSideStds:[1/0,1/0,1/0,1/0],outerSideStdLimit:1/0,outerSideQuads:M,innerSideOk:!1,innerSideStds:[1/0,1/0,1/0,1/0],innerSideStdLimit:1/0,innerSideQuads:v};const L=Y.map(G=>G.mean),O=L.reduce((G,z)=>G+z,0)/L.length,D=B.reduce((G,z)=>G+z.mean,0)/B.length,X=Math.abs(D-O),et=Math.max(...L)-Math.min(...L),V=Y.reduce((G,z)=>G+z.std,0)/Y.length,$=Math.max(0,d),H=Math.max(6,Math.min(20,X*.45)),it=L,at=Y.map(G=>G.std),E=Math.max(H,Math.min(30,H*m)),J=B.map(G=>G.std),Q=J.every(G=>G<=E),T=et<=$&&V<=H,j=1/(1+et/Math.max(1,$)+V/Math.max(1,H));return{ok:T,score:j,meanSpread:et,meanSpreadLimit:$,avgStd:V,avgStdLimit:H,contrast:X,outerMean:O,outerSideMeans:it,outerSideStds:at,outerSideStdLimit:H,outerSideQuads:M,innerSideOk:Q,innerSideStds:J,innerSideStdLimit:E,innerSideQuads:v}}function ge(n,t,e,i,r,s,o){const a=r.x-i.x,l=r.y-i.y,u=Math.hypot(a,l);if(!Number.isFinite(u)||u<=1e-6)return null;const h=a/u,c=l/u,m=-c,d=h,f=(i.x+r.x)*.5,p=(i.y+r.y)*.5,g=xe({p1:i,p2:r},o+2);if(!Mt(Gt(g||[i,r],2),t,e))return null;const x=Math.max(8,Math.round(s*2)+1),b=Math.max(8,Math.round(o*2)+1),_=x>1?s*2/(x-1):0,P=b>1?o*2/(b-1):0,S=Array.from({length:x},()=>new Array(b).fill(0));for(let M=0;M<x;M++){const v=-s+_*M;for(let I=0;I<b;I++){const R=-o+P*I,U=f+v*h+R*m,Y=p+v*c+R*d;S[M][I]=dl(n,t,e,U,Y)}}const w=Ps(S,-o,-s,P,_,!0);if(w.length<8)return null;const C=w.map(M=>{const v=M.x,I=M.y;return{x:f+I*h+v*m,y:p+I*c+v*d,weight:M.weight}}),A=ve(C);if(!A)return null;let k=A.dirX,F=A.dirY;return k*h+F*c<0&&(k=-k,F=-F),{line:{p1:{x:A.pointX-k*s,y:A.pointY-F*s},p2:{x:A.pointX+k*s,y:A.pointY+F*s}},fitPoints:C.map(M=>({x:M.x,y:M.y}))}}function dl(n,t,e,i,r){if(t<=0||e<=0||n.length!==t*e)return 0;const s=Math.max(0,Math.min(t-1,i)),o=Math.max(0,Math.min(e-1,r)),a=Math.floor(s),l=Math.floor(o),u=Math.min(t-1,a+1),h=Math.min(e-1,l+1),c=s-a,m=o-l,d=n[l*t+a],f=n[l*t+u],p=n[h*t+a],g=n[h*t+u],y=d*(1-c)+f*c,x=p*(1-c)+g*c;return y*(1-m)+x*m}function fl(n,t,e,i,r){if(i<=0||r<=0||i>=t-1||r>=e-1)return{gx:0,gy:0};const s=r*t+i;return{gx:(n[s+1]-n[s-1])*.5,gy:(n[s+t]-n[s-t])*.5}}function pl(n){if(n.length<20)return null;const t=n.map(w=>Math.max(0,w.weight));let e=0;for(const w of t)e=Math.max(e,w);if(!(e>0))return null;for(let w=0;w<t.length;w++)t[w]/=e;const i=w=>{let C=0,A=0,k=0;for(let H=0;H<n.length;H++){const it=w[H];it>0&&(C+=it,A+=n[H].x*it,k+=n[H].y*it)}if(!(C>0))return null;A/=C,k/=C;let F=0,M=0,v=0;for(let H=0;H<n.length;H++){const it=w[H];if(!(it>0))continue;const at=n[H].x-A,E=n[H].y-k;F+=it*at*at,M+=it*at*E,v+=it*E*E}F/=C,M/=C,v/=C;const I=F+v,R=F*v-M*M,U=-I,Y=R,B=Math.max(0,U*U-4*Y),L=-.5*(U+(U>=0?1:-1)*Math.sqrt(B)),O=Math.abs(L)>1e-12?L:0,D=Math.abs(L)>1e-12?Y/L:I,X=Math.max(O,D);let et=0,V=1;Math.abs(M)>1e-10?(et=X-v,V=M):F>v&&(et=1,V=0);const $=Math.atan2(-et,V);return{centroid:{x:A,y:k},angle:$,totalWeight:C}},r=i(t);if(!r)return null;const s=Math.cos(r.angle),o=Math.sin(r.angle),a=new Array(2*16*8).fill(0);for(let w=0;w<n.length;w++){const C=n[w].x-r.centroid.x,A=n[w].y-r.centroid.y,k=C*s+A*o,F=Math.round(k*8+16*8);if(F>=3&&F<a.length-3)for(let M=-3;M<=3;M++)a[F+M]+=t[w]}let l=16*8;for(let w=-5*8+16*8;w<=5*8+16*8;w++)a[w]>a[l]&&(l=w);let u=l-1;for(;u>1&&a[u]>.05*a[l];)u--;let h=l+1;for(;h<a.length-1&&a[h]>.05*a[l];)h++;let c=Math.max(1,u-8);for(;c>1&&a[c]<=a[u];)c--;let m=Math.min(a.length-1,h+8);for(;m<a.length-1&&a[m]<=a[h];)m++;const d=a.slice();for(let w=1;w<d.length;w++)d[w]+=d[w-1];const f=d[d.length-1];if(!(f>0))return null;let p=0;for(let w=1;w<d.length;w++)Math.abs(d[w]-.1*f)<Math.abs(d[p]-.1*f)&&(p=w);let g=d.length-1;for(let w=d.length-2;w>0;w--)Math.abs(d[w]-.9*f)<Math.abs(d[g]-.9*f)&&(g=w);let y=p/8-16,x=g/8-16;const b=x-y;y-=b*.7,x+=b*.7,y=Math.max((c+u)/16-16,y),x=Math.min((m+h)/16-16,x);const _=t.slice();for(let w=0;w<n.length;w++){const C=n[w].x-r.centroid.x,A=n[w].y-r.centroid.y,k=C*s+A*o;_[w]=k>=y&&k<=x?t[w]**4*(1/(10+Math.abs(k))):0}const P=i(_);if(!P)return null;const S=[];for(let w=0;w<n.length;w++)_[w]>0&&S.push({x:n[w].x,y:n[w].y,weight:_[w]});return S.length<8?null:{centroid:P.centroid,angle:P.angle,keptSamples:S}}function ml(n,t,e,i,r,s=Lt){var C;const o=r.x-i.x,a=r.y-i.y,l=Math.hypot(o,a);if(!Number.isFinite(l)||l<=12)return null;const u=o/l,h=a/l,c=5,m=4*s+.5,d=(A,k,F,M,v)=>{const I={x:A.x-k*l*.5,y:A.y-F*l*.5},R={p1:I,p2:{x:I.x+k*l,y:I.y+F*l}},U=xe(R,m+2),Y=Mt(Gt(U??[R.p1,R.p2],3),t,e),B=[],L=new Map;if(!Y)return{reduced:null,scanlines:L};for(let O=Y.y;O<Y.y+Y.h;O++)for(let D=Y.x;D<Y.x+Y.w;D++){const X=D,et=O,V=X-I.x,$=et-I.y,H=V*k+$*F;if(!(H>c&&H<l-c))continue;const it=X-A.x,at=et-A.y,E=it*M+at*v;if(Math.abs(E)<12){const{gx:J,gy:Q}=fl(n,t,e,D,O),T=J*J+Q*Q;T>0&&B.push({x:X,y:et,weight:T})}if(Math.abs(E)<m){const J=L.get(O);J?(D<J.start&&(J.start=D),D>J.end&&(J.end=D)):L.set(O,{start:D,end:D})}}return{reduced:pl(B),scanlines:L}};let f={x:(i.x+r.x)*.5,y:(i.y+r.y)*.5},p=u,g=h,y=-g,x=p,b=d(f,p,g,y,x);if(!b.reduced)return null;f=b.reduced.centroid,y=Math.cos(b.reduced.angle),x=Math.sin(b.reduced.angle),p=-x,g=y,p*u+g*h<0&&(p=-p,g=-g,y=-y,x=-x);let _=d(f,p,g,y,x);if(!_.reduced)return null;const P=Math.hypot(_.reduced.centroid.x-f.x,_.reduced.centroid.y-f.y);f=_.reduced.centroid,y=Math.cos(_.reduced.angle),x=Math.sin(_.reduced.angle),p=-x,g=y,p*u+g*h<0&&(p=-p,g=-g,y=-y,x=-x);let S=_;if(P>1){const A=d(f,p,g,y,x);A.reduced&&(S=A,f=A.reduced.centroid,y=Math.cos(A.reduced.angle),x=Math.sin(A.reduced.angle),p=-x,g=y,p*u+g*h<0&&(p=-p,g=-g))}const w=(((C=S.reduced)==null?void 0:C.keptSamples)??[]).map(A=>({x:A.x,y:A.y}));return w.length<8?null:{line:{p1:{x:f.x-p*l*.5,y:f.y-g*l*.5},p2:{x:f.x+p*l*.5,y:f.y+g*l*.5}},fitPoints:w,correctedScanlines:S.scanlines}}function gl(n,t,e){var b;const i=n.length,r=((b=n[0])==null?void 0:b.length)??0;if(r===0||i===0)return 0;const s=Math.max(0,Math.min(r-1,t)),o=Math.max(0,Math.min(i-1,e)),a=Math.floor(s),l=Math.floor(o),u=Math.min(r-1,a+1),h=Math.min(i-1,l+1),c=s-a,m=o-l,d=n[l][a],f=n[l][u],p=n[h][a],g=n[h][u],y=d*(1-c)+f*c,x=p*(1-c)+g*c;return y*(1-m)+x*m}function yl(n,t,e,i,r,s,o){var B;const a=i.p2.x-i.p1.x,l=i.p2.y-i.p1.y,u=Math.hypot(a,l);if(!Number.isFinite(u)||u<=1e-6)return null;const h=a/u,c=l/u,m=-c,d=h,f=(i.p1.x+i.p2.x)*.5,p=(i.p1.y+i.p2.y)*.5,g=xe(i,s+2),y=Mt(Gt(g||[i.p1,i.p2],3),t,e);if(!y)return null;const x=ko(n,t,e,y,0,0,(o==null?void 0:o.bayerPattern)||"RGGB",o==null?void 0:o.greenPhase,o==null?void 0:o.blackLevel),b=x.length;if((((B=x[0])==null?void 0:B.length)??0)<6||b<6)return null;const P=Math.max(8,Math.round(r*2)+1),S=Math.max(8,Math.round(s*2)+1),w=P>1?r*2/(P-1):0,C=S>1?s*2/(S-1):0,A=Array.from({length:P},()=>new Array(S).fill(0));for(let L=0;L<P;L++){const O=-r+w*L;for(let D=0;D<S;D++){const X=-s+C*D,et=f+O*h+X*m,V=p+O*c+X*d;A[L][D]=gl(x,et-y.x,V-y.y)}}const{gx:k,gy:F}=Fo(A),M=k>=F,v=Ps(A,-s,-r,C,w,M);if(v.length<8)return null;const I=v.map(L=>{const O=L.x,D=L.y;return{x:f+D*h+O*m,y:p+D*c+O*d,weight:L.weight}}),R=ve(I);if(!R)return null;let U=R.dirX,Y=R.dirY;return U*h+Y*c<0&&(U=-U,Y=-Y),{line:{p1:{x:R.pointX-U*r,y:R.pointY-Y*r},p2:{x:R.pointX+U*r,y:R.pointY+Y*r}},fitPoints:I.map(L=>({x:L.x,y:L.y}))}}function Ce(n){const t=Math.max(0,Math.min(1,n));return t<=.04045?t/12.92:Math.pow((t+.055)/1.055,2.4)}function ar(n,t,e){if(t<=0||e<=0||n.length!==t*e)return new Uint8Array(Math.max(0,t*e));let i=1/0,r=-1/0;for(let f=0;f<n.length;f++){const p=n[f];Number.isFinite(p)&&(p<i&&(i=p),p>r&&(r=p))}if(!Number.isFinite(i)||!Number.isFinite(r)||r<=i+1e-9)return new Uint8Array(t*e);const s=1024,o=new Uint32Array(s),a=r-i;for(let f=0;f<n.length;f++){const p=Math.max(0,Math.min(1,(n[f]-i)/a)),g=Math.min(s-1,Math.max(0,Math.floor(p*(s-1))));o[g]++}const l=n.length,u=f=>{const p=l*f;let g=0;for(let y=0;y<s;y++)if(g+=o[y],g>=p)return i+y/Math.max(1,s-1)*a;return r},h=u(.01),c=u(.99),m=Math.max(1e-9,c-h),d=new Uint8Array(t*e);for(let f=0;f<n.length;f++){const p=Math.max(0,Math.min(1,(n[f]-h)/m));d[f]=Math.round(p*255)}return d}function Qn(n,t,e=0){const i=new Float32Array(n.width*n.height),r=n.data;for(let s=0,o=0;s<r.length;s+=4,o++)i[o]=Ts(r,s,t,e);return ar(i,n.width,n.height)}function xl(n){return Number.isFinite(n)?Math.max(0,Math.min(65535,Number(n))):0}function Ts(n,t,e,i=0){let r=n[t]/255,s=n[t+1]/255,o=n[t+2]/255;e&&(r=Ce(r),s=Ce(s),o=Ce(o));const a=.2126*r+.7152*s+.0722*o;return Math.max(0,a-xl(i)/65535)}function bl(n,t){const e=n.width,i=n.height,r=n.data;if(r.length<e*i*3)return new Uint8Array(e*i);const s=new Float32Array(e*i);for(let o=0;o<e*i;o++){const a=o*3;t!==void 0?s[o]=r[a+t]:s[o]=.2126*r[a]+.7152*r[a+1]+.0722*r[a+2]}return ar(s,e,i)}function _l(n){const t=new Float32Array(n.width*n.height);for(let e=0;e<n.data.length;e++)t[e]=n.data[e];return ar(t,n.width,n.height)}function Qr(n,t,e){const i=Mt(t,n.width,n.height);if(!i)return null;const r=new Uint16Array(i.w*i.h*3),s=n.data;let o=0;for(let a=i.y;a<i.y+i.h;a++)for(let l=i.x;l<i.x+i.w;l++){const u=(a*n.width+l)*4;let h=s[u]/255,c=s[u+1]/255,m=s[u+2]/255;e&&(h=Ce(h),c=Ce(c),m=Ce(m)),r[o++]=Math.max(0,Math.min(65535,Math.round(h*65535))),r[o++]=Math.max(0,Math.min(65535,Math.round(c*65535))),r[o++]=Math.max(0,Math.min(65535,Math.round(m*65535)))}return{data:r,width:i.w,height:i.h}}function qr(n,t,e,i=0){const r=Mt(t,n.width,n.height);if(!r)return null;const s=new Uint16Array(r.w*r.h*3),o=n.data;let a=0;for(let l=r.y;l<r.y+r.h;l++)for(let u=r.x;u<r.x+r.w;u++){const h=(l*n.width+u)*4,c=Math.max(0,Math.min(65535,Math.round(Ts(o,h,e,i)*65535)));s[a++]=c,s[a++]=c,s[a++]=c}return{data:s,width:r.w,height:r.h}}function Ml(n,t){const e=Mt(t,n.width,n.height);if(!e)return null;const i=new Uint16Array(e.w*e.h),r=n.data;let s=0;for(let o=e.y;o<e.y+e.h;o++)for(let a=e.x;a<e.x+e.w;a++){const l=(o*n.width+a)*4;i[s++]=Math.max(0,Math.min(65535,Math.round((.2126*r[l]+.7152*r[l+1]+.0722*r[l+2])*257)))}return{data:i,width:e.w,height:e.h}}function wl(n,t){const e=Mt(t,n.width,n.height);if(!e)return null;const i=new Uint16Array(e.w*e.h);let r=0;for(let s=e.y;s<e.y+e.h;s++){const o=s*n.width;for(let a=e.x;a<e.x+e.w;a++)i[r++]=n.data[o+a]}return{data:i,width:e.w,height:e.h}}function Qt(n,t,e){return{x:n.x*t,y:n.y*e}}function Sl(n,t,e){return{p1:Qt(n.p1,t,e),p2:Qt(n.p2,t,e)}}function ji(n,t){const e=t(n);return{x:Number.isFinite(e.x)?e.x:n.x,y:Number.isFinite(e.y)?e.y:n.y}}function vl(n,t){return n.map(e=>ji(e,t))}function me(n,t,e,i=0,r=0){if(!n||n.length<8)return;const s=n.map(o=>({x:o.x*t-i,y:o.y*e-r})).filter(o=>Number.isFinite(o.x)&&Number.isFinite(o.y));return s.length>=8?s:void 0}function Cl(n,t){return{p1:ji(n.p1,t),p2:ji(n.p2,t)}}function cn(n,t,e){return{p1:{x:n.p1.x-t,y:n.p1.y-e},p2:{x:n.p2.x-t,y:n.p2.y-e}}}function Pl(n,t,e,i){const r=Math.max(0,Math.min(n.width-1,t)),o=(Math.max(0,Math.min(n.height-1,e))*n.width+r)*4;let a=n.data[o]/255,l=n.data[o+1]/255,u=n.data[o+2]/255;return i&&(a=Ce(a),l=Ce(l),u=Ce(u)),(.2126*a+.7152*l+.0722*u)*65535}function Is(n){return n.kind==="u16-mono"}function en(n){return n.width}function nn(n){return n.height}function ii(n,t,e,i){if(Is(n)){const r=Math.max(0,Math.min(n.width-1,t)),s=Math.max(0,Math.min(n.height-1,e));return n.data[s*n.width+r]}return Pl(n,t,e,i)}function kl(n,t,e,i){if(Is(n)&&n.coordinateSpace==="distorted-padded"){const r=Math.round(n.paddingOffsetX??0),s=Math.round(n.paddingOffsetY??0);return ii(n,t+r,e+s,i)}return ii(n,t,e,i)}function Fl(n,t,e,i,r=3){const o=[...n,{x:(n[0].x+n[1].x+n[2].x+n[3].x)*.25,y:(n[0].y+n[1].y+n[2].y+n[3].y)*.25},{x:(n[0].x+n[1].x)*.5,y:(n[0].y+n[1].y)*.5},{x:(n[1].x+n[2].x)*.5,y:(n[1].y+n[2].y)*.5},{x:(n[2].x+n[3].x)*.5,y:(n[2].y+n[3].y)*.5},{x:(n[3].x+n[0].x)*.5,y:(n[3].y+n[0].y)*.5}].map(a=>Pe(a,t)).filter(a=>Number.isFinite(a.x)&&Number.isFinite(a.y));return o.length===0?null:Mt(Gt(o,r),e,i)}function or(n,t,e,i){const r=new Map;for(let s=n.y;s<n.y+n.h;s++)for(let o=n.x;o<n.x+n.w;o++){const a=Ut({x:o,y:s},t);if(!Number.isFinite(a.x)||!Number.isFinite(a.y))continue;const l=Math.round(a.x),u=Math.round(a.y);if(l<0||u<0||l>=e||u>=i)continue;const h=r.get(u);h?(l<h.start&&(h.start=l),l>h.end&&(h.end=l)):r.set(u,{start:l,end:l})}return r}function Al(n,t,e,i,r,s,o){const a=new Map,l=t.p2.x-t.p1.x,u=t.p2.y-t.p1.y,h=Math.hypot(l,u);if(!Number.isFinite(h)||h<=1e-6)return a;const c=l/h,m=u/h,d=-m,f=c,p={x:(t.p1.x+t.p2.x)*.5,y:(t.p1.y+t.p2.y)*.5},g=Math.max(1,e+1),y=Math.max(1,i+1.5);for(let x=n.y;x<n.y+n.h;x++)for(let b=n.x;b<n.x+n.w;b++){const _=b+.5,P=x+.5,S=_-p.x,w=P-p.y,C=S*c+w*m;if(!Number.isFinite(C)||Math.abs(C)>g)continue;const A=S*d+w*f;if(!Number.isFinite(A)||Math.abs(A)>y)continue;const k=Ut({x:_,y:P},r);if(!Number.isFinite(k.x)||!Number.isFinite(k.y))continue;const F=Math.round(k.x),M=Math.round(k.y);if(F<0||M<0||F>=s||M>=o)continue;const v=a.get(M);v?(F<v.start&&(v.start=F),F>v.end&&(v.end=F)):a.set(M,{start:F,end:F})}return a}function Ns(n,t,e,i,r,s){const o=new Map,a=t.p2.x-t.p1.x,l=t.p2.y-t.p1.y,u=Math.hypot(a,l);if(!Number.isFinite(u)||u<=1e-6)return o;const h=a/u,c=l/u,m=-c,d=h,f={x:(t.p1.x+t.p2.x)*.5,y:(t.p1.y+t.p2.y)*.5},p=Math.max(1,e+1),g=Math.max(1,i+1.5),y=Mt(n,r,s);if(!y)return o;for(let x=y.y;x<y.y+y.h;x++)for(let b=y.x;b<y.x+y.w;b++){const _=b-f.x,P=x-f.y,S=_*h+P*c;if(!Number.isFinite(S)||Math.abs(S)>p)continue;const w=_*m+P*d;if(!Number.isFinite(w)||Math.abs(w)>g)continue;const C=o.get(x);C?(b<C.start&&(C.start=b),b>C.end&&(C.end=b)):o.set(x,{start:b,end:b})}return o}function Rs(n,t,e,i){const r=new Map;for(const[s,o]of n)for(let a=o.start;a<=o.end;a++){const l=Pe({x:a,y:s},t);if(!Number.isFinite(l.x)||!Number.isFinite(l.y))continue;const u=Math.round(l.x),h=Math.round(l.y);if(u<0||h<0||u>=e||h>=i)continue;const c=r.get(h);c?(u<c.start&&(c.start=u),u>c.end&&(c.end=u)):r.set(h,{start:u,end:u})}return r}function lr(n){return Math.abs(n.k1)<1e-4&&Math.abs(n.k2)<1e-4}function Tl(n){return[{x:n.x,y:n.y},{x:n.x+n.w,y:n.y},{x:n.x+n.w,y:n.y+n.h},{x:n.x,y:n.y+n.h}]}function qt(n,t,e,i,r){return Pe({x:i.x+n*t,y:i.y+n*e},r)}function cr(n,t,e,i,r){const o=qt(n,t,e,i,r),a=qt(n+1e-4,t,e,i,r);return{x:(a.x-o.x)/1e-4,y:(a.y-o.y)/1e-4}}function Ls(n,t,e,i,r,s){let o=.01;const a=h=>{const c=qt(h,t,e,i,s);return Math.hypot(c.x-r.x,c.y-r.y)},l=a(n),u=a(n+o);if(!Number.isFinite(l)||!Number.isFinite(u))return null;if(l>u){let h=n,c=n+o;for(let m=0;m<24;m++){o*=2;const d=h+o,f=a(d),p=a(c);if(!Number.isFinite(f)||!Number.isFinite(p))break;if(f>=p)return{a:h,b:d};h=c,c=d}}else{let h=n,c=n+o;for(let m=0;m<24;m++){o*=2;const d=c-o,f=a(d),p=a(h);if(!Number.isFinite(f)||!Number.isFinite(p))break;if(f>=p)return{a:d,b:c};c=h,h=d}}return{a:n-Math.max(.5,o),b:n+Math.max(.5,o)}}function Il(n,t,e=33){const i=n.p2.x-n.p1.x,r=n.p2.y-n.p1.y,s=Math.hypot(i,r);if(!Number.isFinite(s)||s<=1e-6)return[Pe(n.p1,t),Pe(n.p2,t)];const o=i/s,a=r/s,l={x:(n.p1.x+n.p2.x)*.5,y:(n.p1.y+n.p2.y)*.5},u=s*.5,h=Math.max(9,e),c=[];for(let m=0;m<h;m++){const d=h===1?.5:m/(h-1),f=-u+d*(u*2);c.push(qt(f,o,a,l,t))}return c}function Nl(n,t,e,i,r,s,o=1){const a=n.p2.x-n.p1.x,l=n.p2.y-n.p1.y,u=Math.hypot(a,l);if(!Number.isFinite(u)||u<=1e-6)return null;const h=a/u,c=l/u,m={x:(n.p1.x+n.p2.x)*.5,y:(n.p1.y+n.p2.y)*.5},d=Math.max(24,Math.round(e*2)+1),f=[];for(let p=0;p<d;p++){const g=d===1?.5:p/(d-1),y=-e+g*(e*2),x=qt(y,h,c,m,t),b=cr(y,h,c,m,t),_=Math.hypot(b.x,b.y);if(!Number.isFinite(_)||_<=1e-9)continue;const P=-b.y/_,S=b.x/_;f.push({x:x.x+P*(i+o),y:x.y+S*(i+o)},{x:x.x-P*(i+o),y:x.y-S*(i+o)})}if(f.length<2){const p={p1:Pe(n.p1,t),p2:Pe(n.p2,t)},g=xe(p,i+o);return g?Mt(Gt(g,o),r,s):null}return Mt(Gt(f,o),r,s)}function tn(n,t,e,i,r){return Ut({x:i.x+n*t,y:i.y+n*e},r)}function Rl(n,t,e,i,r,s){let o=.01;const a=h=>{const c=tn(h,t,e,i,s);return Math.hypot(c.x-r.x,c.y-r.y)},l=a(n),u=a(n+o);if(!Number.isFinite(l)||!Number.isFinite(u))return null;if(l>u){let h=n,c=n+o;for(let m=0;m<24;m++){o*=2;const d=h+o,f=a(d),p=a(c);if(!Number.isFinite(f)||!Number.isFinite(p))break;if(f>=p)return{a:h,b:d};h=c,c=d}}else{let h=n,c=n+o;for(let m=0;m<24;m++){o*=2;const d=c-o,f=a(d),p=a(h);if(!Number.isFinite(f)||!Number.isFinite(p))break;if(f>=p)return{a:d,b:c};c=h,h=d}}return{a:n-Math.max(.5,o),b:n+Math.max(.5,o)}}function ur(n,t,e){const i=(t.x-n.x)*(t.y-e.y)-(t.x-e.x)*(t.y-n.y);if(Math.abs(i)<=1e-12)return .5*(n.x+e.x);const r=(t.x-n.x)*(t.x-n.x)*(t.y-e.y)-(t.x-e.x)*(t.x-e.x)*(t.y-n.y),s=t.x-.5*r/i;return Number.isFinite(s)?s:.5*(n.x+e.x)}function Ll(n,t,e,i,r){const o=tn(n,t,e,i,r),a=tn(n+1e-4,t,e,i,r);return{x:(a.x-o.x)/1e-4,y:(a.y-o.y)/1e-4}}function El(n,t,e,i,r,s,o=!1,a){const l=[],u=[],h=a?Lt*2:Lt,c=Math.max(1,Math.min(s,h)),m=i.p2.x-i.p1.x,d=i.p2.y-i.p1.y,f=Math.hypot(m,d);if(!Number.isFinite(f)||f<=1e-6)return null;const p=m/f,g=d/f,y={x:(i.p1.x+i.p2.x)*.5,y:(i.p1.y+i.p2.y)*.5},x={p1:Ut(i.p1,e),p2:Ut(i.p2,e)},b=x.p2.x-x.p1.x,_=x.p2.y-x.p1.y,P=Math.hypot(b,_);if(!Number.isFinite(P)||P<=1e-6)return null;const S=b/P,w=_/P,C=-w,A=S,k={x:(x.p1.x+x.p2.x)*.5,y:(x.p1.y+x.p2.y)*.5},F=Mt(a||Gt(xe(i,s+2)??[i.p1,i.p2],2),n.width,n.height);if(!F)return null;const M=Al(F,i,r,c,e,en(t),nn(t));if(M.size===0)return null;const v=!lr(e);for(const[L,O]of M)if(!(L<0||L>=nn(t)))for(let D=O.start;D<=O.end;D++){if(D<0||D>=en(t))continue;const X={x:D,y:L};let et,V;if(v){const $=Pe(X,e);if(!Number.isFinite($.x)||!Number.isFinite($.y)||Math.round($.x)<0||Math.round($.x)>=n.width||Math.round($.y)<0||Math.round($.y)>=n.height)continue;const H=$.x-y.x,it=$.y-y.y,at=H*p+it*g;if(!Number.isFinite(at))continue;et=at,V=H*-g+it*p;const E=Rl(at,p,g,y,X,e);if(!E)continue;const J=.5*(E.a+E.b),Q=tn(E.a,p,g,y,e),T=tn(J,p,g,y,e),j=tn(E.b,p,g,y,e),G=ur({x:E.a,y:Math.hypot(Q.x-X.x,Q.y-X.y)},{x:J,y:Math.hypot(T.x-X.x,T.y-X.y)},{x:E.b,y:Math.hypot(j.x-X.x,j.y-X.y)});if(!Number.isFinite(G))continue;et=G;const z=Ll(G,p,g,y,e),Z=Math.hypot(z.x,z.y);if(!Number.isFinite(Z)||Z<=1e-9)continue;const tt=z.x/Z,ct=-(z.y/Z),dt=tt,rt=tn(G,p,g,y,e);V=(X.x-rt.x)*ct+(X.y-rt.y)*dt}else{const $=X.x-k.x,H=X.y-k.y;et=$*S+H*w,V=$*C+H*A}!Number.isFinite(et)||Math.abs(et)>r||!Number.isFinite(V)||Math.abs(V)>c||(l.push(V),u.push(ii(t,D,L,o)))}if(l.length<8)return null;const I=Ut(i.p1,e),R=Ut(i.p2,e),U=R.x-I.x,Y=R.y-I.y,B=Math.abs(U)>=Math.abs(Y)?1:2;return bn(l,u,B,h)}function Ul(n,t,e,i,r,s,o=!1,a,l,u,h=!1){const c=[],m=[],d=a?Lt*2:Lt,f=Math.max(1,Math.min(s,d)),p=i.p2.x-i.p1.x,g=i.p2.y-i.p1.y,y=Math.hypot(p,g);if(!Number.isFinite(y)||y<=1e-6)return null;const x=p/y,b=g/y,_=-b,P=x,S={x:(i.p1.x+i.p2.x)*.5,y:(i.p1.y+i.p2.y)*.5},w=Mt(a||Gt(xe(i,d*4+2)??[i.p1,i.p2],2),n.width,n.height);if(!w)return null;const C=l??(u?or(Mt(u,en(t),nn(t))??u,e,n.width,n.height):Ns(w,i,Math.max(1,r),f*4+.5,n.width,n.height));if(C.size===0)return null;const A=Rs(C,e,en(t),nn(t));if(A.size===0)return null;const k=!lr(e);for(const[M,v]of A)for(let I=v.start;I<=v.end;I++){const R={x:I,y:M},U=Ut(R,e);if(!Number.isFinite(U.x)||!Number.isFinite(U.y)||Math.round(U.x)<0||Math.round(U.x)>=n.width||Math.round(U.y)<0||Math.round(U.y)>=n.height)continue;const Y=U.x-S.x,B=U.y-S.y,L=Y*x+B*b;let O=Y*_+B*P;if(k){const D=Ls(L,x,b,S,R,e);if(!D)continue;const X=.5*(D.a+D.b),et=qt(D.a,x,b,S,e),V=qt(X,x,b,S,e),$=qt(D.b,x,b,S,e),H=ur({x:D.a,y:Math.hypot(et.x-R.x,et.y-R.y)},{x:X,y:Math.hypot(V.x-R.x,V.y-R.y)},{x:D.b,y:Math.hypot($.x-R.x,$.y-R.y)});if(!Number.isFinite(H))continue;const it=cr(H,x,b,S,e),at=Math.hypot(it.x,it.y);if(!Number.isFinite(at)||at<=1e-9)continue;const E=it.x/at,Q=-(it.y/at),T=E,j=qt(H,x,b,S,e);O=(R.x-j.x)*Q+(R.y-j.y)*T}!Number.isFinite(L)||Math.abs(L)>Math.max(1,r)||!Number.isFinite(O)||Math.abs(O)>f||(c.push(O),m.push(kl(t,I,M,o)))}if(c.length<8)return null;const F=Math.abs(p)>=Math.abs(g)?1:2;return h?An(c,m,F,d):bn(c,m,F,d)}function Bl(n,t,e,i,r,s){const o=n.width,a=n.height,l=(n.bayerPattern||"RGGB").toUpperCase(),u=s!=null&&s.correctedRect?Lt*2:Lt,h=Math.max(1,Math.min(r,u)),c=(s==null?void 0:s.restrictToStrip)??!0,m=e.p2.x-e.p1.x,d=e.p2.y-e.p1.y,f=Math.hypot(m,d);if(!Number.isFinite(f)||f<=1e-6)return null;const p=m/f,g=d/f,y=-g,x=p,b={x:(e.p1.x+e.p2.x)*.5,y:(e.p1.y+e.p2.y)*.5},_={p1:{x:b.x-p*Math.max(1,i),y:b.y-g*Math.max(1,i)},p2:{x:b.x+p*Math.max(1,i),y:b.y+g*Math.max(1,i)}},P=xe(_,h+2),S=(s!=null&&s.fixedRawRect?Mt(s.fixedRawRect,o,a):null)??(s!=null&&s.correctedRect?rr(s.correctedRect,t,o,a):null)??(P?Fl(P,t,o,a,2):null);if(!S)return null;const w=[],C=[];for(let F=S.y;F<S.y+S.h;F++){const M=F*o;for(let v=S.x;v<S.x+S.w;v++){if(!vt(v,F,l,s==null?void 0:s.greenPhase))continue;const I=Ut({x:v,y:F},t);if(!Number.isFinite(I.x)||!Number.isFinite(I.y))continue;const R=I.x-b.x,U=I.y-b.y,Y=R*p+U*g;if(!Number.isFinite(Y)||c&&Math.abs(Y)>Math.max(1,i))continue;const B=R*y+U*x;if(!Number.isFinite(B)||c&&Math.abs(B)>h)continue;w.push(B);let L;L=Math.max(0,n.data[M+v]-xn(s==null?void 0:s.blackLevel,v,F)),C.push(L)}}if(w.length<8)return null;const A=Math.abs(m)>=Math.abs(d)?1:2,k=Math.max(2,(s==null?void 0:s.shortSidePxOverride)??(c?h*2:Math.min(S.w,S.h)));return hi(w,C,k,s==null?void 0:s.manualBinSize,A,s==null?void 0:s.preferAutoPerEdgeBin,!1,!!(s!=null&&s.forceLegacyModel))}function Es(n,t,e,i,r,s=!1,o,a,l,u=!1){if(a&&l){const C=l.p2.x-l.p1.x,A=l.p2.y-l.p1.y,k=Math.hypot(C,A);if(Number.isFinite(k)&&k>1e-6)return El(a,n,t,l,Math.max(1,k*.5),r,s,o)}const h=[],c=[],m=Lt,d=e.p2.x-e.p1.x,f=e.p2.y-e.p1.y,p=Math.hypot(d,f);if(!Number.isFinite(p)||p<=1e-6)return null;const g=d/p,y=f/p,x={x:(e.p1.x+e.p2.x)*.5,y:(e.p1.y+e.p2.y)*.5},b=(o?Mt(o,en(n),nn(n)):null)??Nl(e,t,i+1,r+1,en(n),nn(n),1);if(!b)return null;const _=or(b,t,en(n),nn(n));if(_.size===0)return null;const P=-y,S=g;for(const[C,A]of _)for(let k=A.start;k<=A.end;k++){const F={x:k,y:C};let M=(F.x-x.x)*g+(F.y-x.y)*y,v=(F.x-x.x)*P+(F.y-x.y)*S;!Number.isFinite(M)||Math.abs(M)>i+1||!Number.isFinite(v)||Math.abs(v)>=m||(h.push(v),c.push(ii(n,k,C,s)))}if(h.length<8)return null;const w=Math.abs(d)>=Math.abs(f)?1:2;return u?An(h,c,w,m):bn(h,c,w,m)}function Dl(n,t,e){var m;const i=e.sourceMode??(t.isThreePlane?"three-plane":"rggb-raw"),r=e.useQuadraticProjection!==!1,s=!!e.forceRenderedMeasurement,o=n.width,a=n.height,l=e.threePlaneChannel,u=sr(e.detectionTuning),h=e.monochromeBlackLevel??0;if(i==="rggb-raw"&&!s){if(!t||t.isThreePlane)return null;const d=Ho(t,e.greenPhase),f=o/Math.max(1,t.width),p=a/Math.max(1,t.height);return{sourceMode:i,detectionGray:d,detectionWidth:t.width,detectionHeight:t.height,detectToMeasureX:1,detectToMeasureY:1,detectToDisplayX:f,detectToDisplayY:p,measureToDisplayX:f,measureToDisplayY:p,detectPointToDisplay:g=>Qt(g,f,p),measurePointToDisplay:g=>Qt(g,f,p),displayPointToDetect:g=>Qt(g,1/Math.max(1e-9,f),1/Math.max(1e-9,p)),measureUsesDisplayLine:!1,measureWidth:t.width,measureHeight:t.height,refineLine:(g,y,x)=>yl(t.data,t.width,t.height,{p1:g,p2:y},x*.5,Math.max(4,x*.2),{greenPhase:e.greenPhase,bayerPattern:t.bayerPattern})||ge(d,t.width,t.height,g,y,x*.5,Math.max(4,x*.2)),measureEdge:(g,y,x,b,_)=>Se(t.data,t.width,t.height,g,y,x,b,{greenOnly:!0,greenPhase:e.greenPhase,bayerPattern:t.bayerPattern,blackLevel:e.blackLevel??void 0,preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?me(_==null?void 0:_.fitPoints,1,1):void 0})}}if(s){const d=!!e.distortionCurveApplied&&!!e.distortionModel,f=i==="rggb-raw"&&!!e.distortionCorrected&&!!e.distortionModel&&!t.isThreePlane,p=!!e.distortionCorrected&&!!e.distortionModel&&!!e.distortionOriginalSamplingPlane,g=!!e.distortionCorrected&&!!e.distortionSamplingPlane,y=n,x=Qn(y,!!e.sfrHasGamma,i==="unmix-bw"?h:0);return{sourceMode:i,detectionGray:x,detectionWidth:y.width,detectionHeight:y.height,detectToMeasureX:1,detectToMeasureY:1,detectToDisplayX:1,detectToDisplayY:1,measureToDisplayX:1,measureToDisplayY:1,detectPointToDisplay:b=>b,measurePointToDisplay:b=>b,displayPointToDetect:b=>b,measureUsesDisplayLine:!1,measureWidth:y.width,measureHeight:y.height,refineLine:(b,_,P)=>(p?ml(x,n.width,n.height,b,_,Lt):null)||ge(x,n.width,n.height,b,_,P*.5,Math.max(4,P*.2)),measureEdge:(b,_,P,S,w)=>{const C=e.distortionModel?rr(b,e.distortionModel,t.width,t.height):null;if(f){const M={p1:Ut(_.p1,e.distortionModel),p2:Ut(_.p2,e.distortionModel)},v=Math.hypot(M.p2.x-M.p1.x,M.p2.y-M.p1.y),I=Math.max(2,v*.5*u.sampleHalfWidthRatio);return co(t,e.distortionModel,M,Math.max(1,v*.5),I,{greenPhase:e.greenPhase,blackLevel:e.blackLevel??void 0,correctedRect:b})}if(f)return Bl(t,e.distortionModel,_,P,S,{greenPhase:e.greenPhase,blackLevel:e.blackLevel??void 0,correctedRect:b,fixedRawRect:C,preferAutoPerEdgeBin:!0});if(p)return Ul(n,e.distortionOriginalSamplingPlane,e.distortionModel,_,P,S,!!e.sfrHasGamma,b,(w==null?void 0:w.correctedScanlines)??null,C);if(d){const M={p1:Ut(_.p1,e.distortionModel),p2:Ut(_.p2,e.distortionModel)},v=Math.hypot(M.p2.x-M.p1.x,M.p2.y-M.p1.y);return Es(e.distortionOriginalSamplingPlane??e.distortionSamplingPlane??e.distortionSamplingImage??n,e.distortionModel,M,Math.max(1,v*.5),S,!!e.sfrHasGamma,b,e.distortionBaseImage??n,_)}if(g){const M=wl(e.distortionSamplingPlane,b);if(!M)return null;const v=cn(_,b.x,b.y);return Se(M.data,M.width,M.height,{x:0,y:0,w:M.width,h:M.height},v,P,S,{preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?me(w==null?void 0:w.fitPoints,1,1,b.x,b.y):void 0})}const k=i==="unmix-bw"?qr(n,b,!!e.sfrHasGamma,h):Qr(n,b,!!e.sfrHasGamma);if(!k)return null;const F=cn(_,b.x,b.y);return Se(k.data,k.width,k.height,{x:0,y:0,w:k.width,h:k.height},F,P,S,{isThreePlane:!0,threePlaneChannel:void 0,preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?me(w==null?void 0:w.fitPoints,1,1,b.x,b.y):void 0})}}}if(i==="three-plane"){if(t.isThreePlane&&!e.sfrHasGamma){const f=bl(t,l),p=o/Math.max(1,t.width),g=a/Math.max(1,t.height);return{sourceMode:i,detectionGray:f,detectionWidth:t.width,detectionHeight:t.height,detectToMeasureX:1,detectToMeasureY:1,detectToDisplayX:p,detectToDisplayY:g,measureToDisplayX:p,measureToDisplayY:g,detectPointToDisplay:y=>Qt(y,p,g),measurePointToDisplay:y=>Qt(y,p,g),displayPointToDetect:y=>Qt(y,1/Math.max(1e-9,p),1/Math.max(1e-9,g)),measureUsesDisplayLine:!1,measureWidth:t.width,measureHeight:t.height,refineLine:(y,x,b)=>ge(f,t.width,t.height,y,x,b*.5,Math.max(4,b*.2)),measureEdge:(y,x,b,_,P)=>Se(t.data,t.width,t.height,y,x,b,_,{isThreePlane:!0,threePlaneChannel:l,preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?me(P==null?void 0:P.fitPoints,1,1):void 0})}}const d=Qn(n,!!e.sfrHasGamma);return{sourceMode:i,detectionGray:d,detectionWidth:n.width,detectionHeight:n.height,detectToMeasureX:1,detectToMeasureY:1,detectToDisplayX:1,detectToDisplayY:1,measureToDisplayX:1,measureToDisplayY:1,detectPointToDisplay:f=>f,measurePointToDisplay:f=>f,displayPointToDetect:f=>f,measureUsesDisplayLine:!1,measureWidth:n.width,measureHeight:n.height,refineLine:(f,p,g)=>ge(d,n.width,n.height,f,p,g*.5,Math.max(4,g*.2)),measureEdge:(f,p,g,y,x)=>{const b=Qr(n,f,!!e.sfrHasGamma);if(!b)return null;const _=cn(p,f.x,f.y);return Se(b.data,b.width,b.height,{x:0,y:0,w:b.width,h:b.height},_,g,y,{isThreePlane:!0,threePlaneChannel:l,preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?me(x==null?void 0:x.fitPoints,1,1,f.x,f.y):void 0})}}}if(i==="unmix-bw"){if(t&&!t.isThreePlane&&e.displaySettings){const f=Wa(t,e.displaySettings,e.blackLevel??e.monochromeBlackLevel??void 0);if(f){const p=_l(f),g=o/Math.max(1,t.width),y=a/Math.max(1,t.height);return{sourceMode:i,detectionGray:p,detectionWidth:t.width,detectionHeight:t.height,detectToMeasureX:1,detectToMeasureY:1,detectToDisplayX:g,detectToDisplayY:y,measureToDisplayX:g,measureToDisplayY:y,detectPointToDisplay:x=>Qt(x,g,y),measurePointToDisplay:x=>Qt(x,g,y),displayPointToDetect:x=>Qt(x,1/Math.max(1e-9,g),1/Math.max(1e-9,y)),measureUsesDisplayLine:!1,measureWidth:t.width,measureHeight:t.height,refineLine:(x,b,_)=>ge(p,t.width,t.height,x,b,_*.5,Math.max(4,_*.2)),measureEdge:(x,b,_,P,S)=>Se(f.data,f.width,f.height,x,b,_,P,{preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?me(S==null?void 0:S.fitPoints,1,1):void 0})}}}const d=Qn(n,!!e.sfrHasGamma,h);return{sourceMode:i,detectionGray:d,detectionWidth:n.width,detectionHeight:n.height,detectToMeasureX:1,detectToMeasureY:1,detectToDisplayX:1,detectToDisplayY:1,measureToDisplayX:1,measureToDisplayY:1,detectPointToDisplay:f=>f,measurePointToDisplay:f=>f,displayPointToDetect:f=>f,measureUsesDisplayLine:!1,measureWidth:n.width,measureHeight:n.height,refineLine:(f,p,g)=>ge(d,n.width,n.height,f,p,g*.5,Math.max(4,g*.2)),measureEdge:(f,p,g,y,x)=>{const b=qr(n,f,!!e.sfrHasGamma,h);if(!b)return null;const _=cn(p,f.x,f.y);return Se(b.data,b.width,b.height,{x:0,y:0,w:b.width,h:b.height},_,g,y,{isThreePlane:!0,threePlaneChannel:void 0,preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?me(x==null?void 0:x.fitPoints,1,1,f.x,f.y):void 0})}}}const c=Qn(n,!1);if(t&&!t.isThreePlane&&((m=e.displaySettings)==null?void 0:m.renderMode)==="advanced-zero-dep"&&e.displaySettings.advancedZeroDep){const d=o/Math.max(1,t.width),f=a/Math.max(1,t.height);return{sourceMode:i,detectionGray:c,detectionWidth:n.width,detectionHeight:n.height,detectToMeasureX:t.width/Math.max(1,n.width),detectToMeasureY:t.height/Math.max(1,n.height),detectToDisplayX:1,detectToDisplayY:1,measureToDisplayX:d,measureToDisplayY:f,detectPointToDisplay:p=>p,measurePointToDisplay:p=>Qt(p,d,f),displayPointToDetect:p=>p,measureUsesDisplayLine:!1,measureWidth:t.width,measureHeight:t.height,refineLine:(p,g,y)=>ge(c,n.width,n.height,p,g,y*.5,Math.max(4,y*.2)),measureEdge:(p,g,y,x,b)=>{const _=Ya(t,p,e.displaySettings);if(!_||_.width<8||_.height<8)return null;const P=cn(g,p.x,p.y);return Se(_.data,_.width,_.height,{x:0,y:0,w:_.width,h:_.height},P,y,x,{preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?me(b==null?void 0:b.fitPoints,t.width/Math.max(1,n.width),t.height/Math.max(1,n.height),p.x,p.y):void 0})}}}if(t&&!t.isThreePlane){const d=o/Math.max(1,t.width),f=a/Math.max(1,t.height);return{sourceMode:i,detectionGray:c,detectionWidth:n.width,detectionHeight:n.height,detectToMeasureX:t.width/Math.max(1,n.width),detectToMeasureY:t.height/Math.max(1,n.height),detectToDisplayX:1,detectToDisplayY:1,measureToDisplayX:d,measureToDisplayY:f,detectPointToDisplay:p=>p,measurePointToDisplay:p=>Qt(p,d,f),displayPointToDetect:p=>p,measureUsesDisplayLine:!1,measureWidth:t.width,measureHeight:t.height,refineLine:(p,g,y)=>ge(c,n.width,n.height,p,g,y*.5,Math.max(4,y*.2)),measureEdge:(p,g,y,x,b)=>Se(t.data,t.width,t.height,p,g,y,x,{blackLevel:e.blackLevel??void 0,preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?me(b==null?void 0:b.fitPoints,t.width/Math.max(1,n.width),t.height/Math.max(1,n.height)):void 0})}}return{sourceMode:i,detectionGray:c,detectionWidth:n.width,detectionHeight:n.height,detectToMeasureX:1,detectToMeasureY:1,detectToDisplayX:1,detectToDisplayY:1,measureToDisplayX:1,measureToDisplayY:1,detectPointToDisplay:d=>d,measurePointToDisplay:d=>d,displayPointToDetect:d=>d,measureUsesDisplayLine:!1,measureWidth:n.width,measureHeight:n.height,refineLine:(d,f,p)=>ge(c,n.width,n.height,d,f,p*.5,Math.max(4,p*.2)),measureEdge:(d,f,p,g,y)=>{const x=Ml(n,d);if(!x)return null;const b=cn(f,d.x,d.y);return Do(x.data,x.width,x.height,{x:0,y:0,w:x.width,h:x.height},b,p,g,{blackLevel:e.monochromeBlackLevel??void 0,preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?me(y==null?void 0:y.fitPoints,1,1,d.x,d.y):void 0})}}}function Ol(n,t,e){var u,h,c,m,d,f,p,g;if(!n||!t)return[];(u=e.onProgress)==null||u.call(e,"Preparing source context...",0);const i=Dl(n,t,e);if(!i)return[];(h=e.onProgress)==null||h.call(e,"Preparing source context...",.08);const r=sr(e.detectionTuning),s=Math.min(1e3,Math.max(1,e.maxRegions??1e3)),o=Math.max(4,e.maxEdges??s*4);(c=e.onProgress)==null||c.call(e,"Detecting candidates...",.12);const a=al(i.detectionGray,i.detectionWidth,i.detectionHeight,s,e.detectionTuning,(y,x)=>{var b;(b=e.onProgress)==null||b.call(e,y,.12+.08*Math.max(0,Math.min(1,x)))});(m=e.onProgress)==null||m.call(e,"Detecting candidates...",.2);const l=[];for(let y=0;y<a.length;y++){const x=a.length<=0?1:y/a.length;if((d=e.onProgress)==null||d.call(e,`Measuring edges: region ${y+1}/${a.length}`,.2+.72*Math.min(1,x)),l.length>=o)break;const b=a[y],_=b.corners,P=`auto-region-${y+1}`;for(let S=0;S<4&&((f=e.onProgress)==null||f.call(e,`Measuring edges: region ${y+1}/${a.length}, edge ${S+1}/4`,.2+.72*Math.min(1,(y+S/4)/Math.max(1,a.length))),!(l.length>=o));S=S+1){const w=_[S],C=_[(S+1)%4],A=C.x-w.x,k=C.y-w.y,F=Math.hypot(A,k);if(!Number.isFinite(F)||F<24)continue;const M=.125,v={x:w.x+A*M,y:w.y+k*M},I={x:C.x-A*M,y:C.y-k*M},R=Math.hypot(I.x-v.x,I.y-v.y);if(!Number.isFinite(R)||R<12)continue;const U=i.refineLine(v,I,R),Y=(U!=null&&U.fitPoints?Gr(U.fitPoints):null)||(U==null?void 0:U.line)||{p1:v,p2:I},B=Sl(Y,i.detectToMeasureX,i.detectToMeasureY),L=vl((U==null?void 0:U.fitPoints)??[],i.detectPointToDisplay),O=(L.length>=2?Gr(L):null)||Cl(B,i.measurePointToDisplay),D=i.measureUsesDisplayLine?O:B,X=D.p2.x-D.p1.x,et=D.p2.y-D.p1.y,V=Math.hypot(X,et);if(!Number.isFinite(V)||V<=1e-6)continue;const $=O.p2.x-O.p1.x,H=O.p2.y-O.p1.y,it=Math.hypot($,H);if(!Number.isFinite(it)||it<=1e-6)continue;const at=!!e.distortionCurveApplied&&!!e.distortionModel,E=$/it,J=H/it;let Q=J,T=-E;const j=(O.p1.x+O.p2.x)*.5,G=(O.p1.y+O.p2.y)*.5,z=i.detectPointToDisplay({x:b.centerX,y:b.centerY}),Z=z.x,tt=z.y;(j-Z)*Q+(G-tt)*T<0&&(Q=-Q,T=-T);const nt=V*.5,ct=Math.max(2,V*r.sampleHalfWidthRatio),dt=Math.max(2,it*r.sampleHalfWidthRatio),rt=at?{p1:Ut(O.p1,e.distortionModel),p2:Ut(O.p2,e.distortionModel)}:void 0,Bt=rt?Math.max(1,Math.hypot(rt.p2.x-rt.p1.x,rt.p2.y-rt.p1.y)*.5):nt,Vt=at?O:B,W=at?dt:ct,ut=xe(Vt,W);if(!ut)continue;const ot=Mt(Gt(ut,2),at?n.width:i.measureWidth,at?n.height:i.measureHeight);if(!ot)continue;const pt=e.distortionCorrected&&e.distortionModel&&i.sourceMode==="rggb-raw"?rr(ot,e.distortionModel,t.width,t.height):null,mt=at?Es(e.distortionSamplingPlane??e.distortionSamplingImage??n,e.distortionModel,rt,Bt,W,!!e.sfrHasGamma,ot,e.distortionBaseImage??n,Vt):i.measureEdge(ot,D,nt,W,U);if(!mt||(mt.autoLikeUsed=!0,!Uo(mt,e.useDeshading,0)))continue;const te=e.useNR?-1:12,Dt=zo([mt],te,null,e.useDeshading,0,!0);if(!Dt||Dt.mtf50===null||!Bo(Dt.lsfCropped))continue;const _t=rt?Il(rt,e.distortionModel,Math.max(21,Math.round(it*.5))):mt.quadraticProjectionUsed?Eo(L,O,Math.max(21,Math.round(it*.5))):void 0,le=rt&&_t&&_t.length>=2?Gt(_t,dt+2):null,Kt=le?Tl(le):xe(O,dt);if(!Kt)continue;const ke=e.distortionCorrected?ot:le??Gt(Kt,2);let K={x:j+Q*(dt+12),y:G+T*(dt+12)},Ft=Vr(E,J);if(_t&&_t.length>=3){const Tt=Math.floor(_t.length/2),xt=_t[Math.max(0,Tt-1)],wt=_t[Math.min(_t.length-1,Tt+1)],ee=_t[Tt],$t=wt.x-xt.x,_n=wt.y-xt.y,Be=Math.hypot($t,_n);if(Be>1e-6){const De=_n/Be,Oe=-$t/Be;Ft=Vr($t/Be,_n/Be);const Ge={x:ee.x-Z,y:ee.y-tt},Ve=Ge.x*De+Ge.y*Oe>=0?1:-1;K={x:ee.x+De*Ve*(dt+12),y:ee.y+Oe*Ve*(dt+12)}}}l.push({id:`${P}-edge-${S+1}`,regionId:P,sourceMode:i.sourceMode,edgeIndex:S,label:Dt.mtf50.toFixed(3),mtf50:Dt.mtf50,angle:Ft,orientation:mt.orientation,edgeData:mt,sourceRect:ke,rawSourceRect:(i.sourceMode==="rggb-raw"?pt??ot:pt)??void 0,quad:Kt,line:O,originalLine:B,curveBaseLine:rt,curvePoints:_t,labelPoint:K,ridgePoints:L,outerSideMeans:b.outerSideMeans,outerSideQuads:b.outerSideQuads,distortionCorrected:e.distortionCorrected??!1})}}return(p=e.onProgress)==null||p.call(e,"Finalizing results...",.98),(g=e.onProgress)==null||g.call(e,"Finalizing results...",1),l}const Gl=n=>!n.blackLevels||n.blackLevels.length<4?null:[Number(n.blackLevels[0])||0,Number(n.blackLevels[1])||0,Number(n.blackLevels[2])||0,Number(n.blackLevels[3])||0],Di=(n,t)=>{t instanceof ArrayBuffer&&(n.includes(t)||n.push(t))};self.onmessage=async n=>{var o,a;const{id:t,buffer:e,detect:i,options:r}=n.data,s=performance.now();try{const l=performance.now(),u=await za(e),h=performance.now()-l;let c=0,m=[];if(i&&!u.isXTrans){const f=u.isThreePlane?"three-plane":"rggb-raw",p=performance.now();m=Ol({width:u.width,height:u.height},u,{...r,sourceMode:f,forceRenderedMeasurement:!1,blackLevel:(r==null?void 0:r.blackLevel)??Gl(u),onProgress:(g,y)=>{self.postMessage({id:t,type:"progress",stage:g,progress:y})}}),c=performance.now()-p}const d=[];Di(d,e),Di(d,(o=u.data)==null?void 0:o.buffer),Di(d,(a=u.floatData)==null?void 0:a.buffer),self.postMessage({id:t,type:"result",success:!0,raw:u,rawFileBuffer:e,measurements:m,timings:{decodeMs:h,detectMs:c,totalMs:performance.now()-s}},d)}catch(l){self.postMessage({id:t,type:"result",success:!1,error:(l==null?void 0:l.message)||String(l)})}};
