var Os=Object.defineProperty;var Vs=(n,t,e)=>t in n?Os(n,t,{enumerable:!0,configurable:!0,writable:!0,value:e}):n[t]=e;var bt=(n,t,e)=>Vs(n,typeof t!="symbol"?t+"":t,e);function Gs(n){const[[t,e,i],[r,s,o],[a,l,u]]=n,d=t*(s*u-o*l)-e*(r*u-o*a)+i*(r*l-s*a);if(Math.abs(d)<1e-10)return null;const c=1/d;return[[(s*u-o*l)*c,(i*l-e*u)*c,(e*o-i*s)*c],[(o*a-r*u)*c,(t*u-i*a)*c,(i*r-t*o)*c],[(r*l-s*a)*c,(e*a-t*l)*c,(t*s-e*r)*c]]}function Xs(n,t){const[[e,i,r],[s,o,a],[l,u,d]]=n,[c,m,h]=t;return[e*c+i*m+r*h,s*c+o*m+a*h,l*c+u*m+d*h]}class kn{constructor(t){bt(this,"size");bt(this,"isPowerOfTwo");bt(this,"_real");bt(this,"_imag");bt(this,"_scratch",null);bt(this,"_rev",null);bt(this,"_m",0);bt(this,"_internalFFT",null);bt(this,"_chirpReal",null);bt(this,"_chirpImag",null);bt(this,"_bReal",null);bt(this,"_bImag",null);bt(this,"_hanning",null);bt(this,"_windowSumSq",0);this.size=t,this.isPowerOfTwo=(t&t-1)===0&&t>0,this._real=new Float32Array(t),this._imag=new Float32Array(t),this.isPowerOfTwo?this.initRadix2():this.initBluestein()}initRadix2(){const t=this.size,e=Math.log2(t);this._rev=new Uint32Array(t);for(let i=0;i<t;i++){let r=0,s=i;for(let o=0;o<e;o++)r=r<<1|s&1,s>>>=1;this._rev[i]=r}}initBluestein(){const t=this.size;this._m=Math.pow(2,Math.ceil(Math.log2(2*t-1))),this._internalFFT=new kn(this._m),this._chirpReal=new Float32Array(t),this._chirpImag=new Float32Array(t);for(let r=0;r<t;r++){const s=-Math.PI*(r*r)/t;this._chirpReal[r]=Math.cos(s),this._chirpImag[r]=Math.sin(s)}const e=new Float32Array(this._m),i=new Float32Array(this._m);for(let r=0;r<t;r++)e[r]=this._chirpReal[r],i[r]=-this._chirpImag[r];for(let r=1;r<t;r++)e[this._m-r]=e[r],i[this._m-r]=i[r];this._internalFFT.transform(e,i),this._bReal=new Float32Array(this._internalFFT._real),this._bImag=new Float32Array(this._internalFFT._imag)}initHanning(){if(this._hanning)return;const t=this.size;this._hanning=new Float32Array(t);let e=0;for(let i=0;i<t;i++){const r=.5*(1-Math.cos(2*Math.PI*i/(t-1)));this._hanning[i]=r,e+=r*r}this._windowSumSq=e}transform(t,e){this.isPowerOfTwo?this.transformRadix2(t,e):this.transformBluestein(t,e)}transformRadix2(t,e){const i=this.size,r=this._rev,s=this._real,o=this._imag;if(t===s)for(let a=0;a<i;a++){const l=r[a];if(a<l){const u=s[a],d=o[a];s[a]=s[l],o[a]=o[l],s[l]=u,o[l]=d}}else for(let a=0;a<i;a++){const l=r[a];s[a]=t[l],o[a]=e?e[l]:0}for(let a=2;a<=i;a*=2){const l=a/2,u=-2*Math.PI/a,d=Math.cos(u),c=Math.sin(u);for(let m=0;m<i;m+=a){let h=1,f=0;for(let p=0;p<l;p++){const g=m+p,x=m+p+l,y=h*s[x]-f*o[x],b=h*o[x]+f*s[x],_=s[g],F=o[g];s[g]=_+y,o[g]=F+b,s[x]=_-y,o[x]=F-b;const v=h*d-f*c,w=h*c+f*d;h=v,f=w}}}}transformBluestein(t,e){const i=this.size,r=this._m,s=this._internalFFT,o=s._real,a=s._imag;o.fill(0),a.fill(0);for(let c=0;c<i;c++){const m=t[c],h=e?e[c]:0,f=this._chirpReal[c],p=this._chirpImag[c];o[c]=m*f-h*p,a[c]=m*p+h*f}s.transformRadix2(o,a);for(let c=0;c<r;c++){const m=s._real[c],h=s._imag[c],f=this._bReal[c],p=this._bImag[c];s._real[c]=m*f-h*p,s._imag[c]=m*p+h*f}const l=s._real,u=s._imag;for(let c=0;c<r;c++)u[c]=-u[c];s.transformRadix2(l,u);const d=1/r;for(let c=0;c<i;c++){const m=s._real[c]*d,h=-s._imag[c]*d,f=this._chirpReal[c],p=this._chirpImag[c];this._real[c]=m*f-h*p,this._imag[c]=m*p+h*f}}calculateSpectrum(t,e,i=!1){const r=this.size;let s=0;for(let d=0;d<r;d++)s+=t[d];const o=s/r;this._scratch||(this._scratch=new Float32Array(r));const a=this._scratch;if(i){this.initHanning();const d=this._hanning;for(let c=0;c<r;c++)a[c]=(t[c]-o)*d[c]}else for(let d=0;d<r;d++)a[d]=t[d]-o;this.transform(a);const l=e.length;let u=1/r;i&&this._windowSumSq>0&&(u=1/this._windowSumSq);for(let d=0;d<l;d++){const c=this._real[d],m=this._imag[d];e[d]+=(c*c+m*m)*u}}calculateSpectrumWindow(t,e,i,r=!1){const s=this.size;let o=0;for(let c=0;c<s;c++)o+=t[e+c];const a=o/s;this._scratch||(this._scratch=new Float32Array(s));const l=this._scratch;if(r){this.initHanning();const c=this._hanning;for(let m=0;m<s;m++)l[m]=(t[e+m]-a)*c[m]}else for(let c=0;c<s;c++)l[c]=t[e+c]-a;this.transform(l);const u=i.length;let d=1/s;r&&this._windowSumSq>0&&(d=1/this._windowSumSq);for(let c=0;c<u;c++){const m=this._real[c],h=this._imag[c];i[c]+=(m*m+h*h)*d}}}const zs={"Sony ILCE-7RM5":"0.82 -0.2976 -0.0719 -0.4296 1.2053 0.2532 -0.0429 0.1282 0.5774"};let Pi=null;async function Ys(n){return Pi||(Pi=(async()=>{if(typeof window.loadPyodide!="function")throw new Error("Pyodide missing: window.loadPyodide not found.");const t=await window.loadPyodide();return await t.loadPackage("numpy"),t})()),Pi}var Ws=`#!/usr/bin/env python3
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
`,js=`#!/usr/bin/env python3
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
`,Hs=`#!/usr/bin/env python3
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
`,Qs="/assets/joraw-BAR7UfD7.js",qs="/assets/joraw-CZg3hVRg.wasm";const pr=512,Ks=zs["Sony ILCE-7RM5"].split(/\s+/).map(Number).filter(Number.isFinite);let Ci=null,ki=null;function $s(n){if(n.byteLength<8)return null;const t=n.getUint16(0,!1);return t===18761?!0:t===19789?!1:null}function jr(n,t,e){const i=Math.min(n.length,t+e);let r="";for(let s=t;s<i;s++){const o=n[s];if(o===0)break;r+=String.fromCharCode(o)}return r.trim()}function Qt(n,t,e,i,r){const s=e===1||e===2||e===7?1:e===3||e===8?2:e===4||e===9?4:0;if(!s)return[];const o=s*i,a=o<=4?r:n.getUint32(r,t);if(a<0||a+o>n.byteLength)return[];const l=[];for(let u=0;u<i;u++){const d=a+u*s;e===1||e===2||e===7?l.push(n.getUint8(d)):e===3?l.push(n.getUint16(d,t)):e===8?l.push(n.getInt16(d,t)):e===4?l.push(n.getUint32(d,t)):e===9&&l.push(n.getInt32(d,t))}return l}function mr(n,t,e,i,r,s){if(i!==2||r<=0)return"";const o=r<=4?s:t.getUint32(s,e);return o<0||o>=n.length?"":jr(n,o,r)}function Hr(n){const t=new Uint8Array(n),e=new DataView(n),i=$s(e);if(i===null||e.getUint16(2,i)!==42)return null;const s=c=>e.getUint16(c,i),o=c=>e.getUint32(c,i),a=[o(4)],l=new Set;let u="",d="";for(;a.length;){const c=a.pop();if(l.has(c)||c<=0||c+2>e.byteLength)continue;l.add(c);const m=s(c);if(c+2+m*12+4>e.byteLength)continue;const h=new Map;for(let _=0;_<m;_++){const F=c+2+_*12,v=s(F),w=s(F+2),P=o(F+4),T=F+8;h.set(v,{type:w,count:P,valueOffset:T})}const f=h.get(271),p=h.get(272);f&&!u&&(u=mr(t,e,i,f.type,f.count,f.valueOffset)),p&&!d&&(d=mr(t,e,i,p.type,p.count,p.valueOffset));const g=h.get(330);if(g){const _=Qt(e,i,g.type,g.count,g.valueOffset);for(const F of _)a.push(F)}const x=h.get(259),y=h.get(262);if(x&&y){const _=Qt(e,i,x.type,x.count,x.valueOffset)[0],F=Qt(e,i,y.type,y.count,y.valueOffset)[0];if(_===32766&&F===32803){const v=Qt(e,i,h.get(256).type,h.get(256).count,h.get(256).valueOffset)[0],w=Qt(e,i,h.get(257).type,h.get(257).count,h.get(257).valueOffset)[0],P=Qt(e,i,h.get(258).type,h.get(258).count,h.get(258).valueOffset)[0],T=Qt(e,i,h.get(273).type,h.get(273).count,h.get(273).valueOffset)[0],C=Qt(e,i,h.get(279).type,h.get(279).count,h.get(279).valueOffset)[0],k=h.get(33422)?Qt(e,i,h.get(33422).type,h.get(33422).count,h.get(33422).valueOffset):[0,1,1,2];h.get(29456)&&Qt(e,i,h.get(29456).type,h.get(29456).count,h.get(29456).valueOffset);const M=h.get(50717)?Qt(e,i,h.get(50717).type,h.get(50717).count,h.get(50717).valueOffset)[0]:16383,S=h.get(50719)?Qt(e,i,h.get(50719).type,h.get(50719).count,h.get(50719).valueOffset):[],R=h.get(50720)?Qt(e,i,h.get(50720).type,h.get(50720).count,h.get(50720).valueOffset):[];if(T+pr+16>t.length||T+C>t.length)return null;const A=T+pr,L=jr(t,A,4),X=t[A+8]<<8|t[A+9],D=t[A+10]<<8|t[A+11],E=t[A+12]<<8|t[A+13],U=t[A+14]<<8|t[A+15],O=E>>4&63,z=U>>13,J=U>>10&3,V=D*2,Q=C>=4?(t[T]|t[T+1]<<8|t[T+2]<<16|t[T+3]<<24)>>>0:0,W=X===v&&V===w;let tt=!1;if(Q>=1&&Q<=16&&C>=8+Q*24){const q=new Map,I=new Map;let H=!0;for(let G=0;G<Q;G++){const Y=T+8+G*24,et=e.getUint32(Y+8,!0),nt=e.getUint32(Y+12,!0),it=e.getUint32(Y+16,!0),ct=e.getUint32(Y+20,!0);if(!it||!ct||et+it>v||nt+ct>w){H=!1;break}const dt=q.get(nt);if(dt!==void 0&&dt!==ct){H=!1;break}q.set(nt,ct),I.set(nt,(I.get(nt)||0)+it)}if(H){const G=Array.from(q.keys()).sort((et,nt)=>et-nt);let Y=0;for(const et of G){if(et!==Y||I.get(et)!==v){H=!1;break}Y+=q.get(et)}tt=H&&Y===w}}const at=Q>=1&&Q<=16&&X>0&&V>0&&v%X===0&&w%V===0&&Q===v/X*(w/V);if(L!=="A000"&&L!=="0000"||!W&&!at&&!tt||O!==16||z!==3||J!==3)return null;const B=[1024,1024,1024,1024],Z=k.slice(0,4).map(q=>q===0?"R":q===2?"B":"G").join("")||"RGGB";return{width:v,height:w,bitsPerSample:P,compression:_,photometric:F,blackLevel:B,whiteLevel:Number(M||16383),cfaPattern:Z,defaultCropOrigin:S.length>=2?[Number(S[0]),Number(S[1])]:void 0,defaultCropSize:R.length>=2?[Number(R[0]),Number(R[1])]:void 0,make:u||"SONY",model:d||"ILCE-7M5"}}}const b=o(c+2+m*12);b&&a.push(b)}return null}async function Js(n){return Ci||(Ci=(async()=>{const t=await Ys();return t.__jtrSonyCrawHqDecoderReady||(await t.FS.mkdirTree("/sony_craw_hq"),await t.FS.writeFile("/sony_craw_hq/llvc3_bitstream_probe.py",Ws),await t.FS.writeFile("/sony_craw_hq/llvc3_entropy.py",js),await t.FS.writeFile("/sony_craw_hq/llvc3_math.py",Hs),await t.runPythonAsync(`
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
`),t.__jtrSonyCrawHqDecoderReady=!0),t})()),Ci}function Zs(n){return Hr(n)}async function ta(){return ki||(ki=import(Qs).then(n=>{const t=n.default;if(typeof t!="function")throw new Error("Sony cRAW HQ LibRaw WASM module is missing its initializer");return t({locateFile:(e,i)=>e.endsWith("joraw.wasm")?qs:i+e})})),ki}function ea(n){return n==="ILCE-7M5"?Ks:null}function Qr(n,t,e,i){var l;if(n.length!==t.width*t.height)throw new Error(`Sony cRAW HQ decoded size mismatch: got ${n.length}, expected ${t.width*t.height}`);const r=t.model||"ILCE-7M5",s=r.startsWith("Sony ")?r:`Sony ${r}`,o=ea(r),a={...i||{},make:t.make||(i==null?void 0:i.camera_make)||"SONY",model:r,camera_make:t.make||(i==null?void 0:i.camera_make)||"SONY",camera_model:r,UniqueCameraModel:s,sourceFormat:e==="libraw-wasm"?"Sony cRAW HQ / LLVC3 (LibRaw WASM)":"Sony cRAW HQ / LLVC3",sonyCrawHq:{...t,decodeBackend:e},color_desc:t.cfaPattern,black_level_per_channel:t.blackLevel,white_level:t.whiteLevel,color_matrix:o&&o.length===9?o:void 0,idata:{filters:2492765332,colors:3},color_data:{...(i==null?void 0:i.color_data)||{},black:1024,cblack_rawpy_style:t.blackLevel,dng_levels:{...((l=i==null?void 0:i.color_data)==null?void 0:l.dng_levels)||{},dng_cblack:t.blackLevel,dng_whitelevel:t.whiteLevel}}};return{data:n,width:t.width,height:t.height,bayerPattern:t.cfaPattern,blackLevels:t.blackLevel,whiteLevel:t.whiteLevel,metadata:a,isThreePlane:!1,isXTrans:!1}}async function na(n,t,e){const i=typeof performance<"u"?performance.now():Date.now(),r=await ta(),s=typeof performance<"u"?performance.now():Date.now(),o=r.LibRaw||r.JoRaw;if(!o)throw new Error("Sony cRAW HQ LibRaw WASM class not found");const a=new o;try{const l=new Uint8Array(n);a.open(l,{});const u=typeof performance<"u"?performance.now():Date.now();let d=null;try{d=a.metadata(!0)}catch(x){console.warn("[Sony cRAW HQ] fast WASM metadata read failed",x)}const c=typeof performance<"u"?performance.now():Date.now(),m=a.getRawImage(),h=typeof performance<"u"?performance.now():Date.now();if(!m||!m.data)throw new Error("Sony cRAW HQ LibRaw WASM returned no raw image");const f=m.data instanceof Uint16Array?m.data:new Uint16Array(m.data.buffer,m.data.byteOffset||0,m.data.byteLength/2),p=typeof performance<"u"?performance.now():Date.now();if(m.width!==t.width||m.height!==t.height)throw new Error(`Sony cRAW HQ LibRaw WASM dimensions mismatch: got ${m.width}x${m.height}, expected ${t.width}x${t.height}`);const g=typeof performance<"u"?performance.now():Date.now();return console.info("[Sony cRAW HQ] fast decode timings",{width:t.width,height:t.height,backend:"libraw-wasm",wasmReadyMs:Math.round(s-i),openMs:Math.round(u-s),metadataMs:Math.round(c-u),unpackMs:Math.round(h-c),copyMs:Math.round(p-h),totalMs:Math.round(g-i)}),{rawImageData:Qr(f,t,"libraw-wasm",d),info:t}}finally{typeof a.delete=="function"?a.delete():typeof a.close=="function"&&a.close()}}async function ia(n,t,e){const i=typeof performance<"u"?performance.now():Date.now(),r=await Js(),s=typeof performance<"u"?performance.now():Date.now(),o=new Uint8Array(n),a=await fetch(new URL("/assets/sony_llvc3_static_lut4096_padded_u16-FsVBk-IV.bin",import.meta.url));if(!a.ok)throw new Error(`Failed to load Sony LLVC3 sample LUT: HTTP ${a.status}`);const l=new Uint8Array(await a.arrayBuffer()),u=typeof performance<"u"?performance.now():Date.now();r.globals.set("jtr_sony_arw_bytes",o),r.globals.set("jtr_sony_lut_bytes",l);const d=await r.runPythonAsync("jtr_decode_sony_craw_hq(jtr_sony_arw_bytes.to_py(), jtr_sony_lut_bytes.to_py())"),c=typeof performance<"u"?performance.now():Date.now(),m=d.toJs();typeof d.destroy=="function"&&d.destroy(),r.globals.delete("jtr_sony_arw_bytes"),r.globals.delete("jtr_sony_lut_bytes");const h=new Uint8Array(m.byteLength);h.set(m);const f=new Uint16Array(h.buffer),p=typeof performance<"u"?performance.now():Date.now(),g=typeof performance<"u"?performance.now():Date.now();return console.info("[Sony cRAW HQ] decode timings",{width:t.width,height:t.height,backend:"pyodide",pyodideReadyMs:Math.round(s-i),lutLoadMs:Math.round(u-s),llvc3DecodeMs:Math.round(c-u),copyMs:Math.round(p-c),totalMs:Math.round(g-i)}),{rawImageData:Qr(f,t,"pyodide"),info:t}}async function ra(n,t){const e=Hr(n);if(!e)return null;try{return await na(n,e,t)}catch(i){return console.warn("[Sony cRAW HQ] fast WASM decode failed; falling back to Pyodide",i),ia(n,e)}}async function sa(n,t){return ra(n,t)}var ii=typeof self<"u"?self:global;const Cn=typeof navigator<"u",aa=Cn&&typeof HTMLImageElement>"u",qn=!(typeof global>"u"||typeof process>"u"||!process.versions||!process.versions.node),ri=ii.Buffer,Vn=ii.BigInt,si=!!ri,oa=n=>n;function Kn(n,t=oa){if(qn)try{return typeof require=="function"?Promise.resolve(t(require(n))):import(n).then(t)}catch{console.warn(`Couldn't load ${n}`)}}let Wi=ii.fetch;const la=n=>Wi=n;if(!ii.fetch){const n=Kn("http",i=>i),t=Kn("https",i=>i),e=(i,{headers:r}={})=>new Promise(async(s,o)=>{let{port:a,hostname:l,pathname:u,protocol:d,search:c}=new URL(i);const m={method:"GET",hostname:l,path:encodeURI(u)+c,headers:r};a!==""&&(m.port=Number(a));const h=(d==="https:"?await t:await n).request(m,f=>{if(f.statusCode===301||f.statusCode===302){let p=new URL(f.headers.location,i).toString();return e(p,{headers:r}).then(s).catch(o)}s({status:f.statusCode,arrayBuffer:()=>new Promise(p=>{let g=[];f.on("data",x=>g.push(x)),f.on("end",()=>p(Buffer.concat(g)))})})});h.on("error",o),h.end()});la(e)}function st(n,t,e){return t in n?Object.defineProperty(n,t,{value:e,enumerable:!0,configurable:!0,writable:!0}):n[t]=e,n}const $n=n=>qr(n)?void 0:n,ca=n=>n!==void 0;function qr(n){return n===void 0||(n instanceof Map?n.size===0:Object.values(n).filter(ca).length===0)}function kt(n){let t=new Error(n);throw delete t.stack,t}function Ke(n){return(n=function(t){for(;t.endsWith("\0");)t=t.slice(0,-1);return t}(n).trim())===""?void 0:n}function Di(n){let t=function(e){let i=0;return e.ifd0.enabled&&(i+=1024),e.exif.enabled&&(i+=2048),e.makerNote&&(i+=2048),e.userComment&&(i+=1024),e.gps.enabled&&(i+=512),e.interop.enabled&&(i+=100),e.ifd1.enabled&&(i+=1024),i+2048}(n);return n.jfif.enabled&&(t+=50),n.xmp.enabled&&(t+=2e4),n.iptc.enabled&&(t+=14e3),n.icc.enabled&&(t+=6e3),t}const Bi=n=>String.fromCharCode.apply(null,n),gr=typeof TextDecoder<"u"?new TextDecoder("utf-8"):void 0;function Kr(n){return gr?gr.decode(n):si?Buffer.from(n).toString("utf8"):decodeURIComponent(escape(Bi(n)))}class jt{static from(t,e){return t instanceof this&&t.le===e?t:new jt(t,void 0,void 0,e)}constructor(t,e=0,i,r){if(typeof r=="boolean"&&(this.le=r),Array.isArray(t)&&(t=new Uint8Array(t)),t===0)this.byteOffset=0,this.byteLength=0;else if(t instanceof ArrayBuffer){i===void 0&&(i=t.byteLength-e);let s=new DataView(t,e,i);this._swapDataView(s)}else if(t instanceof Uint8Array||t instanceof DataView||t instanceof jt){i===void 0&&(i=t.byteLength-e),(e+=t.byteOffset)+i>t.byteOffset+t.byteLength&&kt("Creating view outside of available memory in ArrayBuffer");let s=new DataView(t.buffer,e,i);this._swapDataView(s)}else if(typeof t=="number"){let s=new DataView(new ArrayBuffer(t));this._swapDataView(s)}else kt("Invalid input argument for BufferView: "+t)}_swapArrayBuffer(t){this._swapDataView(new DataView(t))}_swapBuffer(t){this._swapDataView(new DataView(t.buffer,t.byteOffset,t.byteLength))}_swapDataView(t){this.dataView=t,this.buffer=t.buffer,this.byteOffset=t.byteOffset,this.byteLength=t.byteLength}_lengthToEnd(t){return this.byteLength-t}set(t,e,i=jt){return t instanceof DataView||t instanceof jt?t=new Uint8Array(t.buffer,t.byteOffset,t.byteLength):t instanceof ArrayBuffer&&(t=new Uint8Array(t)),t instanceof Uint8Array||kt("BufferView.set(): Invalid data argument."),this.toUint8().set(t,e),new i(this,e,t.byteLength)}subarray(t,e){return e=e||this._lengthToEnd(t),new jt(this,t,e)}toUint8(){return new Uint8Array(this.buffer,this.byteOffset,this.byteLength)}getUint8Array(t,e){return new Uint8Array(this.buffer,this.byteOffset+t,e)}getString(t=0,e=this.byteLength){return Kr(this.getUint8Array(t,e))}getLatin1String(t=0,e=this.byteLength){let i=this.getUint8Array(t,e);return Bi(i)}getUnicodeString(t=0,e=this.byteLength){const i=[];for(let r=0;r<e&&t+r<this.byteLength;r+=2)i.push(this.getUint16(t+r));return Bi(i)}getInt8(t){return this.dataView.getInt8(t)}getUint8(t){return this.dataView.getUint8(t)}getInt16(t,e=this.le){return this.dataView.getInt16(t,e)}getInt32(t,e=this.le){return this.dataView.getInt32(t,e)}getUint16(t,e=this.le){return this.dataView.getUint16(t,e)}getUint32(t,e=this.le){return this.dataView.getUint32(t,e)}getFloat32(t,e=this.le){return this.dataView.getFloat32(t,e)}getFloat64(t,e=this.le){return this.dataView.getFloat64(t,e)}getFloat(t,e=this.le){return this.dataView.getFloat32(t,e)}getDouble(t,e=this.le){return this.dataView.getFloat64(t,e)}getUintBytes(t,e,i){switch(e){case 1:return this.getUint8(t,i);case 2:return this.getUint16(t,i);case 4:return this.getUint32(t,i);case 8:return this.getUint64&&this.getUint64(t,i)}}getUint(t,e,i){switch(e){case 8:return this.getUint8(t,i);case 16:return this.getUint16(t,i);case 32:return this.getUint32(t,i);case 64:return this.getUint64&&this.getUint64(t,i)}}toString(t){return this.dataView.toString(t,this.constructor.name)}ensureChunk(){}}function Oi(n,t){kt(`${n} '${t}' was not loaded, try using full build of exifr.`)}class ji extends Map{constructor(t){super(),this.kind=t}get(t,e){return this.has(t)||Oi(this.kind,t),e&&(t in e||function(i,r){kt(`Unknown ${i} '${r}'.`)}(this.kind,t),e[t].enabled||Oi(this.kind,t)),super.get(t)}keyList(){return Array.from(this.keys())}}var xe=new ji("file parser"),Ct=new ji("segment parser"),_e=new ji("file reader");function ua(n,t){return typeof n=="string"?yr(n,t):Cn&&!aa&&n instanceof HTMLImageElement?yr(n.src,t):n instanceof Uint8Array||n instanceof ArrayBuffer||n instanceof DataView?new jt(n):Cn&&n instanceof Blob?Vi(n,t,"blob",sn):void kt("Invalid input argument")}function yr(n,t){return(e=n).startsWith("data:")||e.length>1e4?Gi(n,t,"base64"):qn&&n.includes("://")?Vi(n,t,"url",rn):qn?Gi(n,t,"fs"):Cn?Vi(n,t,"url",rn):void kt("Invalid input argument");var e}async function Vi(n,t,e,i){return _e.has(e)?Gi(n,t,e):i?async function(r,s){let o=await s(r);return new jt(o)}(n,i):void kt(`Parser ${e} is not loaded`)}async function Gi(n,t,e){let i=new(_e.get(e))(n,t);return await i.read(),i}const rn=n=>Wi(n).then(t=>t.arrayBuffer()),sn=n=>new Promise((t,e)=>{let i=new FileReader;i.onloadend=()=>t(i.result||new ArrayBuffer),i.onerror=e,i.readAsArrayBuffer(n)});class ha extends Map{get tagKeys(){return this.allKeys||(this.allKeys=Array.from(this.keys())),this.allKeys}get tagValues(){return this.allValues||(this.allValues=Array.from(this.values())),this.allValues}}function Pt(n,t,e){let i=new ha;for(let[r,s]of e)i.set(r,s);if(Array.isArray(t))for(let r of t)n.set(r,i);else n.set(t,i);return i}function an(n,t,e){let i,r=n.get(t);for(i of e)r.set(i[0],i[1])}const Rt=new Map,le=new Map,Ue=new Map,Ie=["chunked","firstChunkSize","firstChunkSizeNode","firstChunkSizeBrowser","chunkSize","chunkLimit"],mn=["jfif","xmp","icc","iptc","ihdr"],on=["tiff",...mn],St=["ifd0","ifd1","exif","gps","interop"],Re=[...on,...St],Ne=["makerNote","userComment"],gn=["translateKeys","translateValues","reviveValues","multiSegment"],Le=[...gn,"sanitize","mergeOutput","silentErrors"];class $r{get translate(){return this.translateKeys||this.translateValues||this.reviveValues}}class wn extends $r{get needed(){return this.enabled||this.deps.size>0}constructor(t,e,i,r){if(super(),st(this,"enabled",!1),st(this,"skip",new Set),st(this,"pick",new Set),st(this,"deps",new Set),st(this,"translateKeys",!1),st(this,"translateValues",!1),st(this,"reviveValues",!1),this.key=t,this.enabled=e,this.parse=this.enabled,this.applyInheritables(r),this.canBeFiltered=St.includes(t),this.canBeFiltered&&(this.dict=Rt.get(t)),i!==void 0)if(Array.isArray(i))this.parse=this.enabled=!0,this.canBeFiltered&&i.length>0&&this.translateTagSet(i,this.pick);else if(typeof i=="object"){if(this.enabled=!0,this.parse=i.parse!==!1,this.canBeFiltered){let{pick:s,skip:o}=i;s&&s.length>0&&this.translateTagSet(s,this.pick),o&&o.length>0&&this.translateTagSet(o,this.skip)}this.applyInheritables(i)}else i===!0||i===!1?this.parse=this.enabled=i:kt(`Invalid options argument: ${i}`)}applyInheritables(t){let e,i;for(e of gn)i=t[e],i!==void 0&&(this[e]=i)}translateTagSet(t,e){if(this.dict){let i,r,{tagKeys:s,tagValues:o}=this.dict;for(i of t)typeof i=="string"?(r=o.indexOf(i),r===-1&&(r=s.indexOf(Number(i))),r!==-1&&e.add(Number(s[r]))):e.add(i)}else for(let i of t)e.add(i)}finalizeFilters(){!this.enabled&&this.deps.size>0?(this.enabled=!0,Jn(this.pick,this.deps)):this.enabled&&this.pick.size>0&&Jn(this.pick,this.deps)}}var Ot={jfif:!1,tiff:!0,xmp:!1,icc:!1,iptc:!1,ifd0:!0,ifd1:!1,exif:!0,gps:!0,interop:!1,ihdr:void 0,makerNote:!1,userComment:!1,multiSegment:!1,skip:[],pick:[],translateKeys:!0,translateValues:!0,reviveValues:!0,sanitize:!0,mergeOutput:!0,silentErrors:!0,chunked:!0,firstChunkSize:void 0,firstChunkSizeNode:512,firstChunkSizeBrowser:65536,chunkSize:65536,chunkLimit:5},xr=new Map;class ln extends $r{static useCached(t){let e=xr.get(t);return e!==void 0||(e=new this(t),xr.set(t,e)),e}constructor(t){super(),t===!0?this.setupFromTrue():t===void 0?this.setupFromUndefined():Array.isArray(t)?this.setupFromArray(t):typeof t=="object"?this.setupFromObject(t):kt(`Invalid options argument ${t}`),this.firstChunkSize===void 0&&(this.firstChunkSize=Cn?this.firstChunkSizeBrowser:this.firstChunkSizeNode),this.mergeOutput&&(this.ifd1.enabled=!1),this.filterNestedSegmentTags(),this.traverseTiffDependencyTree(),this.checkLoadedPlugins()}setupFromUndefined(){let t;for(t of Ie)this[t]=Ot[t];for(t of Le)this[t]=Ot[t];for(t of Ne)this[t]=Ot[t];for(t of Re)this[t]=new wn(t,Ot[t],void 0,this)}setupFromTrue(){let t;for(t of Ie)this[t]=Ot[t];for(t of Le)this[t]=Ot[t];for(t of Ne)this[t]=!0;for(t of Re)this[t]=new wn(t,!0,void 0,this)}setupFromArray(t){let e;for(e of Ie)this[e]=Ot[e];for(e of Le)this[e]=Ot[e];for(e of Ne)this[e]=Ot[e];for(e of Re)this[e]=new wn(e,!1,void 0,this);this.setupGlobalFilters(t,void 0,St)}setupFromObject(t){let e;for(e of(St.ifd0=St.ifd0||St.image,St.ifd1=St.ifd1||St.thumbnail,Object.assign(this,t),Ie))this[e]=Fi(t[e],Ot[e]);for(e of Le)this[e]=Fi(t[e],Ot[e]);for(e of Ne)this[e]=Fi(t[e],Ot[e]);for(e of on)this[e]=new wn(e,Ot[e],t[e],this);for(e of St)this[e]=new wn(e,Ot[e],t[e],this.tiff);this.setupGlobalFilters(t.pick,t.skip,St,Re),t.tiff===!0?this.batchEnableWithBool(St,!0):t.tiff===!1?this.batchEnableWithUserValue(St,t):Array.isArray(t.tiff)?this.setupGlobalFilters(t.tiff,void 0,St):typeof t.tiff=="object"&&this.setupGlobalFilters(t.tiff.pick,t.tiff.skip,St)}batchEnableWithBool(t,e){for(let i of t)this[i].enabled=e}batchEnableWithUserValue(t,e){for(let i of t){let r=e[i];this[i].enabled=r!==!1&&r!==void 0}}setupGlobalFilters(t,e,i,r=i){if(t&&t.length){for(let o of r)this[o].enabled=!1;let s=br(t,i);for(let[o,a]of s)Jn(this[o].pick,a),this[o].enabled=!0}else if(e&&e.length){let s=br(e,i);for(let[o,a]of s)Jn(this[o].skip,a)}}filterNestedSegmentTags(){let{ifd0:t,exif:e,xmp:i,iptc:r,icc:s}=this;this.makerNote?e.deps.add(37500):e.skip.add(37500),this.userComment?e.deps.add(37510):e.skip.add(37510),i.enabled||t.skip.add(700),r.enabled||t.skip.add(33723),s.enabled||t.skip.add(34675)}traverseTiffDependencyTree(){let{ifd0:t,exif:e,gps:i,interop:r}=this;r.needed&&(e.deps.add(40965),t.deps.add(40965)),e.needed&&t.deps.add(34665),i.needed&&t.deps.add(34853),this.tiff.enabled=St.some(s=>this[s].enabled===!0)||this.makerNote||this.userComment;for(let s of St)this[s].finalizeFilters()}get onlyTiff(){return!mn.map(t=>this[t].enabled).some(t=>t===!0)&&this.tiff.enabled}checkLoadedPlugins(){for(let t of on)this[t].enabled&&!Ct.has(t)&&Oi("segment parser",t)}}function br(n,t){let e,i,r,s,o=[];for(r of t){for(s of(e=Rt.get(r),i=[],e))(n.includes(s[0])||n.includes(s[1]))&&i.push(s[0]);i.length&&o.push([r,i])}return o}function Fi(n,t){return n!==void 0?n:t!==void 0?t:void 0}function Jn(n,t){for(let e of t)n.add(e)}st(ln,"default",Ot);class De{constructor(t){st(this,"parsers",{}),st(this,"output",{}),st(this,"errors",[]),st(this,"pushToErrors",e=>this.errors.push(e)),this.options=ln.useCached(t)}async read(t){this.file=await ua(t,this.options)}setup(){if(this.fileParser)return;let{file:t}=this,e=t.getUint16(0);for(let[i,r]of xe)if(r.canHandle(t,e))return this.fileParser=new r(this.options,this.file,this.parsers),t[i]=!0;this.file.close&&this.file.close(),kt("Unknown file format")}async parse(){let{output:t,errors:e}=this;return this.setup(),this.options.silentErrors?(await this.executeParsers().catch(this.pushToErrors),e.push(...this.fileParser.errors)):await this.executeParsers(),this.file.close&&this.file.close(),this.options.silentErrors&&e.length>0&&(t.errors=e),$n(t)}async executeParsers(){let{output:t}=this;await this.fileParser.parse();let e=Object.values(this.parsers).map(async i=>{let r=await i.parse();i.assignToOutput(t,r)});this.options.silentErrors&&(e=e.map(i=>i.catch(this.pushToErrors))),await Promise.all(e)}async extractThumbnail(){this.setup();let{options:t,file:e}=this,i=Ct.get("tiff",t);var r;if(e.tiff?r={start:0,type:"tiff"}:e.jpeg&&(r=await this.fileParser.getOrFindSegment("tiff")),r===void 0)return;let s=await this.fileParser.ensureSegmentChunk(r),o=this.parsers.tiff=new i(s,t,e),a=await o.extractThumbnail();return e.close&&e.close(),a}}async function ai(n,t){let e=new De(t);return await e.read(n),e.parse()}var da=Object.freeze({__proto__:null,parse:ai,Exifr:De,fileParsers:xe,segmentParsers:Ct,fileReaders:_e,tagKeys:Rt,tagValues:le,tagRevivers:Ue,createDictionary:Pt,extendDictionary:an,fetchUrlAsArrayBuffer:rn,readBlobAsArrayBuffer:sn,chunkedProps:Ie,otherSegments:mn,segments:on,tiffBlocks:St,segmentsAndBlocks:Re,tiffExtractables:Ne,inheritables:gn,allFormatters:Le,Options:ln});class oi{constructor(t,e,i){st(this,"errors",[]),st(this,"ensureSegmentChunk",async r=>{let s=r.start,o=r.size||65536;if(this.file.chunked)if(this.file.available(s,o))r.chunk=this.file.subarray(s,o);else try{r.chunk=await this.file.readChunk(s,o)}catch(a){kt(`Couldn't read segment: ${JSON.stringify(r)}. ${a.message}`)}else this.file.byteLength>s+o?r.chunk=this.file.subarray(s,o):r.size===void 0?r.chunk=this.file.subarray(s):kt("Segment unreachable: "+JSON.stringify(r));return r.chunk}),this.extendOptions&&this.extendOptions(t),this.options=t,this.file=e,this.parsers=i}injectSegment(t,e){this.options[t].enabled&&this.createParser(t,e)}createParser(t,e){let i=new(Ct.get(t))(e,this.options,this.file);return this.parsers[t]=i}createParsers(t){for(let e of t){let{type:i,chunk:r}=e,s=this.options[i];if(s&&s.enabled){let o=this.parsers[i];o&&o.append||o||this.createParser(i,r)}}}async readSegments(t){let e=t.map(this.ensureSegmentChunk);await Promise.all(e)}}class oe{static findPosition(t,e){let i=t.getUint16(e+2)+2,r=typeof this.headerLength=="function"?this.headerLength(t,e,i):this.headerLength,s=e+r,o=i-r;return{offset:e,length:i,headerLength:r,start:s,size:o,end:s+o}}static parse(t,e={}){return new this(t,new ln({[this.type]:e}),t).parse()}normalizeInput(t){return t instanceof jt?t:new jt(t)}constructor(t,e={},i){st(this,"errors",[]),st(this,"raw",new Map),st(this,"handleError",r=>{if(!this.options.silentErrors)throw r;this.errors.push(r.message)}),this.chunk=this.normalizeInput(t),this.file=i,this.type=this.constructor.type,this.globalOptions=this.options=e,this.localOptions=e[this.type],this.canTranslate=this.localOptions&&this.localOptions.translate}translate(){this.canTranslate&&(this.translated=this.translateBlock(this.raw,this.type))}get output(){return this.translated?this.translated:this.raw?Object.fromEntries(this.raw):void 0}translateBlock(t,e){let i=Ue.get(e),r=le.get(e),s=Rt.get(e),o=this.options[e],a=o.reviveValues&&!!i,l=o.translateValues&&!!r,u=o.translateKeys&&!!s,d={};for(let[c,m]of t)a&&i.has(c)?m=i.get(c)(m):l&&r.has(c)&&(m=this.translateValue(m,r.get(c))),u&&s.has(c)&&(c=s.get(c)||c),d[c]=m;return d}translateValue(t,e){return e[t]||e.DEFAULT||t}assignToOutput(t,e){this.assignObjectToOutput(t,this.constructor.type,e)}assignObjectToOutput(t,e,i){if(this.globalOptions.mergeOutput)return Object.assign(t,i);t[e]?Object.assign(t[e],i):t[e]=i}}st(oe,"headerLength",4),st(oe,"type",void 0),st(oe,"multiSegment",!1),st(oe,"canHandle",()=>!1);function fa(n){return n===192||n===194||n===196||n===219||n===221||n===218||n===254}function pa(n){return n>=224&&n<=239}function ma(n,t,e){for(let[i,r]of Ct)if(r.canHandle(n,t,e))return i}class _r extends oi{constructor(...t){super(...t),st(this,"appSegments",[]),st(this,"jpegSegments",[]),st(this,"unknownSegments",[])}static canHandle(t,e){return e===65496}async parse(){await this.findAppSegments(),await this.readSegments(this.appSegments),this.mergeMultiSegments(),this.createParsers(this.mergedAppSegments||this.appSegments)}setupSegmentFinderArgs(t){t===!0?(this.findAll=!0,this.wanted=new Set(Ct.keyList())):(t=t===void 0?Ct.keyList().filter(e=>this.options[e].enabled):t.filter(e=>this.options[e].enabled&&Ct.has(e)),this.findAll=!1,this.remaining=new Set(t),this.wanted=new Set(t)),this.unfinishedMultiSegment=!1}async findAppSegments(t=0,e){this.setupSegmentFinderArgs(e);let{file:i,findAll:r,wanted:s,remaining:o}=this;if(!r&&this.file.chunked&&(r=Array.from(s).some(a=>{let l=Ct.get(a),u=this.options[a];return l.multiSegment&&u.multiSegment}),r&&await this.file.readWhole()),t=this.findAppSegmentsInRange(t,i.byteLength),!this.options.onlyTiff&&i.chunked){let a=!1;for(;o.size>0&&!a&&(i.canReadNextChunk||this.unfinishedMultiSegment);){let{nextChunkOffset:l}=i,u=this.appSegments.some(d=>!this.file.available(d.offset||d.start,d.length||d.size));if(a=t>l&&!u?!await i.readNextChunk(t):!await i.readNextChunk(l),(t=this.findAppSegmentsInRange(t,i.byteLength))===void 0)return}}}findAppSegmentsInRange(t,e){e-=2;let i,r,s,o,a,l,{file:u,findAll:d,wanted:c,remaining:m,options:h}=this;for(;t<e;t++)if(u.getUint8(t)===255){if(i=u.getUint8(t+1),pa(i)){if(r=u.getUint16(t+2),s=ma(u,t,r),s&&c.has(s)&&(o=Ct.get(s),a=o.findPosition(u,t),l=h[s],a.type=s,this.appSegments.push(a),!d&&(o.multiSegment&&l.multiSegment?(this.unfinishedMultiSegment=a.chunkNumber<a.chunkCount,this.unfinishedMultiSegment||m.delete(s)):m.delete(s),m.size===0)))break;h.recordUnknownSegments&&(a=oe.findPosition(u,t),a.marker=i,this.unknownSegments.push(a)),t+=r+1}else if(fa(i)){if(r=u.getUint16(t+2),i===218&&h.stopAfterSos!==!1)return;h.recordJpegSegments&&this.jpegSegments.push({offset:t,length:r,marker:i}),t+=r+1}}return t}mergeMultiSegments(){if(!this.appSegments.some(e=>e.multiSegment))return;let t=function(e,i){let r,s,o,a=new Map;for(let l=0;l<e.length;l++)r=e[l],s=r[i],a.has(s)?o=a.get(s):a.set(s,o=[]),o.push(r);return Array.from(a)}(this.appSegments,"type");this.mergedAppSegments=t.map(([e,i])=>{let r=Ct.get(e,this.options);return r.handleMultiSegments?{type:e,chunk:r.handleMultiSegments(i)}:i[0]})}getSegment(t){return this.appSegments.find(e=>e.type===t)}async getOrFindSegment(t){let e=this.getSegment(t);return e===void 0&&(await this.findAppSegments(0,[t]),e=this.getSegment(t)),e}}st(_r,"type","jpeg"),xe.set("jpeg",_r);const ga=[void 0,1,1,2,4,8,1,1,2,4,8,4,8,4];class ya extends oe{parseHeader(){var t=this.chunk.getUint16();t===18761?this.le=!0:t===19789&&(this.le=!1),this.chunk.le=this.le,this.headerParsed=!0}parseTags(t,e,i=new Map){let{pick:r,skip:s}=this.options[e];r=new Set(r);let o=r.size>0,a=s.size===0,l=this.chunk.getUint16(t);t+=2;for(let u=0;u<l;u++){let d=this.chunk.getUint16(t);if(o){if(r.has(d)&&(i.set(d,this.parseTag(t,d,e)),r.delete(d),r.size===0))break}else!a&&s.has(d)||i.set(d,this.parseTag(t,d,e));t+=12}return i}parseTag(t,e,i){let{chunk:r}=this,s=r.getUint16(t+2),o=r.getUint32(t+4),a=ga[s];if(a*o<=4?t+=8:t=r.getUint32(t+8),(s<1||s>13)&&kt(`Invalid TIFF value type. block: ${i.toUpperCase()}, tag: ${e.toString(16)}, type: ${s}, offset ${t}`),t>r.byteLength&&kt(`Invalid TIFF value offset. block: ${i.toUpperCase()}, tag: ${e.toString(16)}, type: ${s}, offset ${t} is outside of chunk size ${r.byteLength}`),s===1)return r.getUint8Array(t,o);if(s===2)return Ke(r.getString(t,o));if(s===7)return r.getUint8Array(t,o);if(o===1)return this.parseTagValue(s,t);{let l=new(function(d){switch(d){case 1:return Uint8Array;case 3:return Uint16Array;case 4:return Uint32Array;case 5:return Array;case 6:return Int8Array;case 8:return Int16Array;case 9:return Int32Array;case 10:return Array;case 11:return Float32Array;case 12:return Float64Array;default:return Array}}(s))(o),u=a;for(let d=0;d<o;d++)l[d]=this.parseTagValue(s,t),t+=u;return l}}parseTagValue(t,e){let{chunk:i}=this;switch(t){case 1:return i.getUint8(e);case 3:return i.getUint16(e);case 4:return i.getUint32(e);case 5:return i.getUint32(e)/i.getUint32(e+4);case 6:return i.getInt8(e);case 8:return i.getInt16(e);case 9:return i.getInt32(e);case 10:return i.getInt32(e)/i.getInt32(e+4);case 11:return i.getFloat(e);case 12:return i.getDouble(e);case 13:return i.getUint32(e);default:kt(`Invalid tiff type ${t}`)}}}class Ai extends ya{static canHandle(t,e){return t.getUint8(e+1)===225&&t.getUint32(e+4)===1165519206&&t.getUint16(e+8)===0}async parse(){this.parseHeader();let{options:t}=this;return t.ifd0.enabled&&await this.parseIfd0Block(),t.exif.enabled&&await this.safeParse("parseExifBlock"),t.gps.enabled&&await this.safeParse("parseGpsBlock"),t.interop.enabled&&await this.safeParse("parseInteropBlock"),t.ifd1.enabled&&await this.safeParse("parseThumbnailBlock"),this.createOutput()}safeParse(t){let e=this[t]();return e.catch!==void 0&&(e=e.catch(this.handleError)),e}findIfd0Offset(){this.ifd0Offset===void 0&&(this.ifd0Offset=this.chunk.getUint32(4))}findIfd1Offset(){if(this.ifd1Offset===void 0){this.findIfd0Offset();let t=this.chunk.getUint16(this.ifd0Offset),e=this.ifd0Offset+2+12*t;this.ifd1Offset=this.chunk.getUint32(e)}}parseBlock(t,e){let i=new Map;return this[e]=i,this.parseTags(t,e,i),i}async parseIfd0Block(){if(this.ifd0)return;let{file:t}=this;this.findIfd0Offset(),this.ifd0Offset<8&&kt("Malformed EXIF data"),!t.chunked&&this.ifd0Offset>t.byteLength&&kt(`IFD0 offset points to outside of file.
this.ifd0Offset: ${this.ifd0Offset}, file.byteLength: ${t.byteLength}`),t.tiff&&await t.ensureChunk(this.ifd0Offset,Di(this.options));let e=this.parseBlock(this.ifd0Offset,"ifd0");return e.size!==0?(this.exifOffset=e.get(34665),this.interopOffset=e.get(40965),this.gpsOffset=e.get(34853),this.xmp=e.get(700),this.iptc=e.get(33723),this.icc=e.get(34675),this.options.sanitize&&(e.delete(34665),e.delete(40965),e.delete(34853),e.delete(700),e.delete(33723),e.delete(34675)),e):void 0}async parseExifBlock(){if(this.exif||(this.ifd0||await this.parseIfd0Block(),this.exifOffset===void 0))return;this.file.tiff&&await this.file.ensureChunk(this.exifOffset,Di(this.options));let t=this.parseBlock(this.exifOffset,"exif");return this.interopOffset||(this.interopOffset=t.get(40965)),this.makerNote=t.get(37500),this.userComment=t.get(37510),this.options.sanitize&&(t.delete(40965),t.delete(37500),t.delete(37510)),this.unpack(t,41728),this.unpack(t,41729),t}unpack(t,e){let i=t.get(e);i&&i.length===1&&t.set(e,i[0])}async parseGpsBlock(){if(this.gps||(this.ifd0||await this.parseIfd0Block(),this.gpsOffset===void 0))return;let t=this.parseBlock(this.gpsOffset,"gps");return t&&t.has(2)&&t.has(4)&&(t.set("latitude",Mr(...t.get(2),t.get(1))),t.set("longitude",Mr(...t.get(4),t.get(3)))),t}async parseInteropBlock(){if(!this.interop&&(this.ifd0||await this.parseIfd0Block(),this.interopOffset!==void 0||this.exif||await this.parseExifBlock(),this.interopOffset!==void 0))return this.parseBlock(this.interopOffset,"interop")}async parseThumbnailBlock(t=!1){if(!this.ifd1&&!this.ifd1Parsed&&(!this.options.mergeOutput||t))return this.findIfd1Offset(),this.ifd1Offset>0&&(this.parseBlock(this.ifd1Offset,"ifd1"),this.ifd1Parsed=!0),this.ifd1}async extractThumbnail(){if(this.headerParsed||this.parseHeader(),this.ifd1Parsed||await this.parseThumbnailBlock(!0),this.ifd1===void 0)return;let t=this.ifd1.get(513),e=this.ifd1.get(514);return this.chunk.getUint8Array(t,e)}get image(){return this.ifd0}get thumbnail(){return this.ifd1}createOutput(){let t,e,i,r={};for(e of St)if(t=this[e],!qr(t))if(i=this.canTranslate?this.translateBlock(t,e):Object.fromEntries(t),this.options.mergeOutput){if(e==="ifd1")continue;Object.assign(r,i)}else r[e]=i;return this.makerNote&&(r.makerNote=this.makerNote),this.userComment&&(r.userComment=this.userComment),r}assignToOutput(t,e){if(this.globalOptions.mergeOutput)Object.assign(t,e);else for(let[i,r]of Object.entries(e))this.assignObjectToOutput(t,i,r)}}function Mr(n,t,e,i){var r=n+t/60+e/3600;return i!=="S"&&i!=="W"||(r*=-1),r}st(Ai,"type","tiff"),st(Ai,"headerLength",10),Ct.set("tiff",Ai);var xa=Object.freeze({__proto__:null,default:da,Exifr:De,fileParsers:xe,segmentParsers:Ct,fileReaders:_e,tagKeys:Rt,tagValues:le,tagRevivers:Ue,createDictionary:Pt,extendDictionary:an,fetchUrlAsArrayBuffer:rn,readBlobAsArrayBuffer:sn,chunkedProps:Ie,otherSegments:mn,segments:on,tiffBlocks:St,segmentsAndBlocks:Re,tiffExtractables:Ne,inheritables:gn,allFormatters:Le,Options:ln,parse:ai});const Hi={ifd0:!1,ifd1:!1,exif:!1,gps:!1,interop:!1,sanitize:!1,reviveValues:!0,translateKeys:!1,translateValues:!1,mergeOutput:!1},Qi=Object.assign({},Hi,{firstChunkSize:4e4,gps:[1,2,3,4]});async function Jr(n){let t=new De(Qi);await t.read(n);let e=await t.parse();if(e&&e.gps){let{latitude:i,longitude:r}=e.gps;return{latitude:i,longitude:r}}}const qi=Object.assign({},Hi,{tiff:!1,ifd1:!0,mergeOutput:!1});async function Zr(n){let t=new De(qi);await t.read(n);let e=await t.extractThumbnail();return e&&si?ri.from(e):e}async function ts(n){let t=await this.thumbnail(n);if(t!==void 0){let e=new Blob([t]);return URL.createObjectURL(e)}}const Ki=Object.assign({},Hi,{firstChunkSize:4e4,ifd0:[274]});async function $i(n){let t=new De(Ki);await t.read(n);let e=await t.parse();if(e&&e.ifd0)return e.ifd0[274]}const Ji=Object.freeze({1:{dimensionSwapped:!1,scaleX:1,scaleY:1,deg:0,rad:0},2:{dimensionSwapped:!1,scaleX:-1,scaleY:1,deg:0,rad:0},3:{dimensionSwapped:!1,scaleX:1,scaleY:1,deg:180,rad:180*Math.PI/180},4:{dimensionSwapped:!1,scaleX:-1,scaleY:1,deg:180,rad:180*Math.PI/180},5:{dimensionSwapped:!0,scaleX:1,scaleY:-1,deg:90,rad:90*Math.PI/180},6:{dimensionSwapped:!0,scaleX:1,scaleY:1,deg:90,rad:90*Math.PI/180},7:{dimensionSwapped:!0,scaleX:1,scaleY:-1,deg:270,rad:270*Math.PI/180},8:{dimensionSwapped:!0,scaleX:1,scaleY:1,deg:270,rad:270*Math.PI/180}});let $e=!0,Je=!0;if(typeof navigator=="object"){let n=navigator.userAgent;if(n.includes("iPad")||n.includes("iPhone")){let t=n.match(/OS (\d+)_(\d+)/);if(t){let[,e,i]=t;$e=Number(e)+.1*Number(i)<13.4,Je=!1}}else if(n.includes("OS X 10")){let[,t]=n.match(/OS X 10[_.](\d+)/);$e=Je=Number(t)<15}if(n.includes("Chrome/")){let[,t]=n.match(/Chrome\/(\d+)/);$e=Je=Number(t)<81}else if(n.includes("Firefox/")){let[,t]=n.match(/Firefox\/(\d+)/);$e=Je=Number(t)<77}}async function es(n){let t=await $i(n);return Object.assign({canvas:$e,css:Je},Ji[t])}class ba extends jt{constructor(...t){super(...t),st(this,"ranges",new _a),this.byteLength!==0&&this.ranges.add(0,this.byteLength)}_tryExtend(t,e,i){if(t===0&&this.byteLength===0&&i){let r=new DataView(i.buffer||i,i.byteOffset,i.byteLength);this._swapDataView(r)}else{let r=t+e;if(r>this.byteLength){let{dataView:s}=this._extend(r);this._swapDataView(s)}}}_extend(t){let e;e=si?ri.allocUnsafe(t):new Uint8Array(t);let i=new DataView(e.buffer,e.byteOffset,e.byteLength);return e.set(new Uint8Array(this.buffer,this.byteOffset,this.byteLength),0),{uintView:e,dataView:i}}subarray(t,e,i=!1){return e=e||this._lengthToEnd(t),i&&this._tryExtend(t,e),this.ranges.add(t,e),super.subarray(t,e)}set(t,e,i=!1){i&&this._tryExtend(e,t.byteLength,t);let r=super.set(t,e);return this.ranges.add(e,r.byteLength),r}async ensureChunk(t,e){this.chunked&&(this.ranges.available(t,e)||await this.readChunk(t,e))}available(t,e){return this.ranges.available(t,e)}}class _a{constructor(){st(this,"list",[])}get length(){return this.list.length}add(t,e,i=0){let r=t+e,s=this.list.filter(o=>wr(t,o.offset,r)||wr(t,o.end,r));if(s.length>0){t=Math.min(t,...s.map(a=>a.offset)),r=Math.max(r,...s.map(a=>a.end)),e=r-t;let o=s.shift();o.offset=t,o.length=e,o.end=r,this.list=this.list.filter(a=>!s.includes(a))}else this.list.push({offset:t,length:e,end:r})}available(t,e){let i=t+e;return this.list.some(r=>r.offset<=t&&i<=r.end)}}function wr(n,t,e){return n<=t&&t<=e}class li extends ba{constructor(t,e){super(0),st(this,"chunksRead",0),this.input=t,this.options=e}async readWhole(){this.chunked=!1,await this.readChunk(this.nextChunkOffset)}async readChunked(){this.chunked=!0,await this.readChunk(0,this.options.firstChunkSize)}async readNextChunk(t=this.nextChunkOffset){if(this.fullyRead)return this.chunksRead++,!1;let e=this.options.chunkSize,i=await this.readChunk(t,e);return!!i&&i.byteLength===e}async readChunk(t,e){if(this.chunksRead++,(e=this.safeWrapAddress(t,e))!==0)return this._readChunk(t,e)}safeWrapAddress(t,e){return this.size!==void 0&&t+e>this.size?Math.max(0,this.size-t):e}get nextChunkOffset(){if(this.ranges.list.length!==0)return this.ranges.list[0].length}get canReadNextChunk(){return this.chunksRead<this.options.chunkLimit}get fullyRead(){return this.size!==void 0&&this.nextChunkOffset===this.size}read(){return this.options.chunked?this.readChunked():this.readWhole()}close(){}}_e.set("blob",class extends li{async readWhole(){this.chunked=!1;let n=await sn(this.input);this._swapArrayBuffer(n)}readChunked(){return this.chunked=!0,this.size=this.input.size,super.readChunked()}async _readChunk(n,t){let e=t?n+t:void 0,i=this.input.slice(n,e),r=await sn(i);return this.set(r,n,!0)}});var Ma=Object.freeze({__proto__:null,default:xa,Exifr:De,fileParsers:xe,segmentParsers:Ct,fileReaders:_e,tagKeys:Rt,tagValues:le,tagRevivers:Ue,createDictionary:Pt,extendDictionary:an,fetchUrlAsArrayBuffer:rn,readBlobAsArrayBuffer:sn,chunkedProps:Ie,otherSegments:mn,segments:on,tiffBlocks:St,segmentsAndBlocks:Re,tiffExtractables:Ne,inheritables:gn,allFormatters:Le,Options:ln,parse:ai,gpsOnlyOptions:Qi,gps:Jr,thumbnailOnlyOptions:qi,thumbnail:Zr,thumbnailUrl:ts,orientationOnlyOptions:Ki,orientation:$i,rotations:Ji,get rotateCanvas(){return $e},get rotateCss(){return Je},rotation:es});_e.set("url",class extends li{async readWhole(){this.chunked=!1;let n=await rn(this.input);n instanceof ArrayBuffer?this._swapArrayBuffer(n):n instanceof Uint8Array&&this._swapBuffer(n)}async _readChunk(n,t){let e=t?n+t-1:void 0,i=this.options.httpHeaders||{};(n||e)&&(i.range=`bytes=${[n,e].join("-")}`);let r=await Wi(this.input,{headers:i}),s=await r.arrayBuffer(),o=s.byteLength;if(r.status!==416)return o!==t&&(this.size=n+o),this.set(s,n,!0)}});jt.prototype.getUint64=function(n){let t=this.getUint32(n),e=this.getUint32(n+4);return t<1048575?t<<32|e:typeof Vn!==void 0?(console.warn("Using BigInt because of type 64uint but JS can only handle 53b numbers."),Vn(t)<<Vn(32)|Vn(e)):void kt("Trying to read 64b value but JS can only handle 53b numbers.")};class wa extends oi{parseBoxes(t=0){let e=[];for(;t<this.file.byteLength-4;){let i=this.parseBoxHead(t);if(e.push(i),i.length===0)break;t+=i.length}return e}parseSubBoxes(t){t.boxes=this.parseBoxes(t.start)}findBox(t,e){return t.boxes===void 0&&this.parseSubBoxes(t),t.boxes.find(i=>i.kind===e)}parseBoxHead(t){let e=this.file.getUint32(t),i=this.file.getString(t+4,4),r=t+8;return e===1&&(e=this.file.getUint64(t+8),r+=8),{offset:t,length:e,kind:i,start:r}}parseBoxFullHead(t){if(t.version!==void 0)return;let e=this.file.getUint32(t.start);t.version=e>>24,t.start+=4}}class ns extends wa{static canHandle(t,e){if(e!==0)return!1;let i=t.getUint16(2);if(i>50)return!1;let r=16,s=[];for(;r<i;)s.push(t.getString(r,4)),r+=4;return s.includes(this.type)}async parse(){let t=this.file.getUint32(0),e=this.parseBoxHead(t);for(;e.kind!=="meta";)t+=e.length,await this.file.ensureChunk(t,16),e=this.parseBoxHead(t);await this.file.ensureChunk(e.offset,e.length),this.parseBoxFullHead(e),this.parseSubBoxes(e),this.options.icc.enabled&&await this.findIcc(e),this.options.tiff.enabled&&await this.findExif(e)}async registerSegment(t,e,i){await this.file.ensureChunk(e,i);let r=this.file.subarray(e,i);this.createParser(t,r)}async findIcc(t){let e=this.findBox(t,"iprp");if(e===void 0)return;let i=this.findBox(e,"ipco");if(i===void 0)return;let r=this.findBox(i,"colr");r!==void 0&&await this.registerSegment("icc",r.offset+12,r.length)}async findExif(t){let e=this.findBox(t,"iinf");if(e===void 0)return;let i=this.findBox(t,"iloc");if(i===void 0)return;let r=this.findExifLocIdInIinf(e),s=this.findExtentInIloc(i,r);if(s===void 0)return;let[o,a]=s;await this.file.ensureChunk(o,a);let l=4+this.file.getUint32(o);o+=l,a-=l,await this.registerSegment("tiff",o,a)}findExifLocIdInIinf(t){this.parseBoxFullHead(t);let e,i,r,s,o=t.start,a=this.file.getUint16(o);for(o+=2;a--;){if(e=this.parseBoxHead(o),this.parseBoxFullHead(e),i=e.start,e.version>=2&&(r=e.version===3?4:2,s=this.file.getString(i+r+2,4),s==="Exif"))return this.file.getUintBytes(i,r);o+=e.length}}get8bits(t){let e=this.file.getUint8(t);return[e>>4,15&e]}findExtentInIloc(t,e){this.parseBoxFullHead(t);let i=t.start,[r,s]=this.get8bits(i++),[o,a]=this.get8bits(i++),l=t.version===2?4:2,u=t.version===1||t.version===2?2:0,d=a+r+s,c=t.version===2?4:2,m=this.file.getUintBytes(i,c);for(i+=c;m--;){let h=this.file.getUintBytes(i,l);i+=l+u+2+o;let f=this.file.getUint16(i);if(i+=2,h===e)return f>1&&console.warn(`ILOC box has more than one extent but we're only processing one
Please create an issue at https://github.com/MikeKovarik/exifr with this file`),[this.file.getUintBytes(i+a,r),this.file.getUintBytes(i+a+r,s)];i+=f*d}}}class is extends ns{}st(is,"type","heic");class Sr extends ns{}st(Sr,"type","avif"),xe.set("heic",is),xe.set("avif",Sr),Pt(Rt,["ifd0","ifd1"],[[256,"ImageWidth"],[257,"ImageHeight"],[258,"BitsPerSample"],[259,"Compression"],[262,"PhotometricInterpretation"],[270,"ImageDescription"],[271,"Make"],[272,"Model"],[273,"StripOffsets"],[274,"Orientation"],[277,"SamplesPerPixel"],[278,"RowsPerStrip"],[279,"StripByteCounts"],[282,"XResolution"],[283,"YResolution"],[284,"PlanarConfiguration"],[296,"ResolutionUnit"],[301,"TransferFunction"],[305,"Software"],[306,"ModifyDate"],[315,"Artist"],[316,"HostComputer"],[317,"Predictor"],[318,"WhitePoint"],[319,"PrimaryChromaticities"],[513,"ThumbnailOffset"],[514,"ThumbnailLength"],[529,"YCbCrCoefficients"],[530,"YCbCrSubSampling"],[531,"YCbCrPositioning"],[532,"ReferenceBlackWhite"],[700,"ApplicationNotes"],[33432,"Copyright"],[33723,"IPTC"],[34665,"ExifIFD"],[34675,"ICC"],[34853,"GpsIFD"],[330,"SubIFD"],[40965,"InteropIFD"],[40091,"XPTitle"],[40092,"XPComment"],[40093,"XPAuthor"],[40094,"XPKeywords"],[40095,"XPSubject"]]),Pt(Rt,"exif",[[33434,"ExposureTime"],[33437,"FNumber"],[34850,"ExposureProgram"],[34852,"SpectralSensitivity"],[34855,"ISO"],[34858,"TimeZoneOffset"],[34859,"SelfTimerMode"],[34864,"SensitivityType"],[34865,"StandardOutputSensitivity"],[34866,"RecommendedExposureIndex"],[34867,"ISOSpeed"],[34868,"ISOSpeedLatitudeyyy"],[34869,"ISOSpeedLatitudezzz"],[36864,"ExifVersion"],[36867,"DateTimeOriginal"],[36868,"CreateDate"],[36873,"GooglePlusUploadCode"],[36880,"OffsetTime"],[36881,"OffsetTimeOriginal"],[36882,"OffsetTimeDigitized"],[37121,"ComponentsConfiguration"],[37122,"CompressedBitsPerPixel"],[37377,"ShutterSpeedValue"],[37378,"ApertureValue"],[37379,"BrightnessValue"],[37380,"ExposureCompensation"],[37381,"MaxApertureValue"],[37382,"SubjectDistance"],[37383,"MeteringMode"],[37384,"LightSource"],[37385,"Flash"],[37386,"FocalLength"],[37393,"ImageNumber"],[37394,"SecurityClassification"],[37395,"ImageHistory"],[37396,"SubjectArea"],[37500,"MakerNote"],[37510,"UserComment"],[37520,"SubSecTime"],[37521,"SubSecTimeOriginal"],[37522,"SubSecTimeDigitized"],[37888,"AmbientTemperature"],[37889,"Humidity"],[37890,"Pressure"],[37891,"WaterDepth"],[37892,"Acceleration"],[37893,"CameraElevationAngle"],[40960,"FlashpixVersion"],[40961,"ColorSpace"],[40962,"ExifImageWidth"],[40963,"ExifImageHeight"],[40964,"RelatedSoundFile"],[41483,"FlashEnergy"],[41486,"FocalPlaneXResolution"],[41487,"FocalPlaneYResolution"],[41488,"FocalPlaneResolutionUnit"],[41492,"SubjectLocation"],[41493,"ExposureIndex"],[41495,"SensingMethod"],[41728,"FileSource"],[41729,"SceneType"],[41730,"CFAPattern"],[41985,"CustomRendered"],[41986,"ExposureMode"],[41987,"WhiteBalance"],[41988,"DigitalZoomRatio"],[41989,"FocalLengthIn35mmFormat"],[41990,"SceneCaptureType"],[41991,"GainControl"],[41992,"Contrast"],[41993,"Saturation"],[41994,"Sharpness"],[41996,"SubjectDistanceRange"],[42016,"ImageUniqueID"],[42032,"OwnerName"],[42033,"SerialNumber"],[42034,"LensInfo"],[42035,"LensMake"],[42036,"LensModel"],[42037,"LensSerialNumber"],[42080,"CompositeImage"],[42081,"CompositeImageCount"],[42082,"CompositeImageExposureTimes"],[42240,"Gamma"],[59932,"Padding"],[59933,"OffsetSchema"],[65e3,"OwnerName"],[65001,"SerialNumber"],[65002,"Lens"],[65100,"RawFile"],[65101,"Converter"],[65102,"WhiteBalance"],[65105,"Exposure"],[65106,"Shadows"],[65107,"Brightness"],[65108,"Contrast"],[65109,"Saturation"],[65110,"Sharpness"],[65111,"Smoothness"],[65112,"MoireFilter"],[40965,"InteropIFD"]]),Pt(Rt,"gps",[[0,"GPSVersionID"],[1,"GPSLatitudeRef"],[2,"GPSLatitude"],[3,"GPSLongitudeRef"],[4,"GPSLongitude"],[5,"GPSAltitudeRef"],[6,"GPSAltitude"],[7,"GPSTimeStamp"],[8,"GPSSatellites"],[9,"GPSStatus"],[10,"GPSMeasureMode"],[11,"GPSDOP"],[12,"GPSSpeedRef"],[13,"GPSSpeed"],[14,"GPSTrackRef"],[15,"GPSTrack"],[16,"GPSImgDirectionRef"],[17,"GPSImgDirection"],[18,"GPSMapDatum"],[19,"GPSDestLatitudeRef"],[20,"GPSDestLatitude"],[21,"GPSDestLongitudeRef"],[22,"GPSDestLongitude"],[23,"GPSDestBearingRef"],[24,"GPSDestBearing"],[25,"GPSDestDistanceRef"],[26,"GPSDestDistance"],[27,"GPSProcessingMethod"],[28,"GPSAreaInformation"],[29,"GPSDateStamp"],[30,"GPSDifferential"],[31,"GPSHPositioningError"]]),Pt(le,["ifd0","ifd1"],[[274,{1:"Horizontal (normal)",2:"Mirror horizontal",3:"Rotate 180",4:"Mirror vertical",5:"Mirror horizontal and rotate 270 CW",6:"Rotate 90 CW",7:"Mirror horizontal and rotate 90 CW",8:"Rotate 270 CW"}],[296,{1:"None",2:"inches",3:"cm"}]]);let Pn=Pt(le,"exif",[[34850,{0:"Not defined",1:"Manual",2:"Normal program",3:"Aperture priority",4:"Shutter priority",5:"Creative program",6:"Action program",7:"Portrait mode",8:"Landscape mode"}],[37121,{0:"-",1:"Y",2:"Cb",3:"Cr",4:"R",5:"G",6:"B"}],[37383,{0:"Unknown",1:"Average",2:"CenterWeightedAverage",3:"Spot",4:"MultiSpot",5:"Pattern",6:"Partial",255:"Other"}],[37384,{0:"Unknown",1:"Daylight",2:"Fluorescent",3:"Tungsten (incandescent light)",4:"Flash",9:"Fine weather",10:"Cloudy weather",11:"Shade",12:"Daylight fluorescent (D 5700 - 7100K)",13:"Day white fluorescent (N 4600 - 5400K)",14:"Cool white fluorescent (W 3900 - 4500K)",15:"White fluorescent (WW 3200 - 3700K)",17:"Standard light A",18:"Standard light B",19:"Standard light C",20:"D55",21:"D65",22:"D75",23:"D50",24:"ISO studio tungsten",255:"Other"}],[37385,{0:"Flash did not fire",1:"Flash fired",5:"Strobe return light not detected",7:"Strobe return light detected",9:"Flash fired, compulsory flash mode",13:"Flash fired, compulsory flash mode, return light not detected",15:"Flash fired, compulsory flash mode, return light detected",16:"Flash did not fire, compulsory flash mode",24:"Flash did not fire, auto mode",25:"Flash fired, auto mode",29:"Flash fired, auto mode, return light not detected",31:"Flash fired, auto mode, return light detected",32:"No flash function",65:"Flash fired, red-eye reduction mode",69:"Flash fired, red-eye reduction mode, return light not detected",71:"Flash fired, red-eye reduction mode, return light detected",73:"Flash fired, compulsory flash mode, red-eye reduction mode",77:"Flash fired, compulsory flash mode, red-eye reduction mode, return light not detected",79:"Flash fired, compulsory flash mode, red-eye reduction mode, return light detected",89:"Flash fired, auto mode, red-eye reduction mode",93:"Flash fired, auto mode, return light not detected, red-eye reduction mode",95:"Flash fired, auto mode, return light detected, red-eye reduction mode"}],[41495,{1:"Not defined",2:"One-chip color area sensor",3:"Two-chip color area sensor",4:"Three-chip color area sensor",5:"Color sequential area sensor",7:"Trilinear sensor",8:"Color sequential linear sensor"}],[41728,{1:"Film Scanner",2:"Reflection Print Scanner",3:"Digital Camera"}],[41729,{1:"Directly photographed"}],[41985,{0:"Normal",1:"Custom",2:"HDR (no original saved)",3:"HDR (original saved)",4:"Original (for HDR)",6:"Panorama",7:"Portrait HDR",8:"Portrait"}],[41986,{0:"Auto",1:"Manual",2:"Auto bracket"}],[41987,{0:"Auto",1:"Manual"}],[41990,{0:"Standard",1:"Landscape",2:"Portrait",3:"Night",4:"Other"}],[41991,{0:"None",1:"Low gain up",2:"High gain up",3:"Low gain down",4:"High gain down"}],[41996,{0:"Unknown",1:"Macro",2:"Close",3:"Distant"}],[42080,{0:"Unknown",1:"Not a Composite Image",2:"General Composite Image",3:"Composite Image Captured While Shooting"}]]);const vr={1:"No absolute unit of measurement",2:"Inch",3:"Centimeter"};Pn.set(37392,vr),Pn.set(41488,vr);const Ti={0:"Normal",1:"Low",2:"High"};function Pr(n){return typeof n=="object"&&n.length!==void 0?n[0]:n}function Cr(n){let t=Array.from(n).slice(1);return t[1]>15&&(t=t.map(e=>String.fromCharCode(e))),t[2]!=="0"&&t[2]!==0||t.pop(),t.join(".")}function Ii(n){if(typeof n=="string"){var[t,e,i,r,s,o]=n.trim().split(/[-: ]/g).map(Number),a=new Date(t,e-1,i);return Number.isNaN(r)||Number.isNaN(s)||Number.isNaN(o)||(a.setHours(r),a.setMinutes(s),a.setSeconds(o)),Number.isNaN(+a)?n:a}}function Sn(n){if(typeof n=="string")return n;let t=[];if(n[1]===0&&n[n.length-1]===0)for(let e=0;e<n.length;e+=2)t.push(kr(n[e+1],n[e]));else for(let e=0;e<n.length;e+=2)t.push(kr(n[e],n[e+1]));return Ke(String.fromCodePoint(...t))}function kr(n,t){return n<<8|t}Pn.set(41992,Ti),Pn.set(41993,Ti),Pn.set(41994,Ti),Pt(Ue,["ifd0","ifd1"],[[50827,function(n){return typeof n!="string"?Kr(n):n}],[306,Ii],[40091,Sn],[40092,Sn],[40093,Sn],[40094,Sn],[40095,Sn]]),Pt(Ue,"exif",[[40960,Cr],[36864,Cr],[36867,Ii],[36868,Ii],[40962,Pr],[40963,Pr]]),Pt(Ue,"gps",[[0,n=>Array.from(n).join(".")],[7,n=>Array.from(n).join(":")]]);class Ri extends oe{static canHandle(t,e){return t.getUint8(e+1)===225&&t.getUint32(e+4)===1752462448&&t.getString(e+4,20)==="http://ns.adobe.com/"}static headerLength(t,e){return t.getString(e+4,34)==="http://ns.adobe.com/xmp/extension/"?79:33}static findPosition(t,e){let i=super.findPosition(t,e);return i.multiSegment=i.extended=i.headerLength===79,i.multiSegment?(i.chunkCount=t.getUint8(e+72),i.chunkNumber=t.getUint8(e+76),t.getUint8(e+77)!==0&&i.chunkNumber++):(i.chunkCount=1/0,i.chunkNumber=-1),i}static handleMultiSegments(t){return t.map(e=>e.chunk.getString()).join("")}normalizeInput(t){return typeof t=="string"?t:jt.from(t).getString()}parse(t=this.chunk){if(!this.localOptions.parse)return t;t=function(s){let o={},a={};for(let l of os)o[l]=[],a[l]=0;return s.replace(Ca,(l,u,d)=>{if(u==="<"){let c=++a[d];return o[d].push(c),`${l}#${c}`}return`${l}#${o[d].pop()}`})}(t);let e=hn.findAll(t,"rdf","Description");e.length===0&&e.push(new hn("rdf","Description",void 0,t));let i,r={};for(let s of e)for(let o of s.properties)i=Pa(o.ns,r),rs(o,i);return function(s){let o;for(let a in s)o=s[a]=$n(s[a]),o===void 0&&delete s[a];return $n(s)}(r)}assignToOutput(t,e){if(this.localOptions.parse)for(let[i,r]of Object.entries(e))switch(i){case"tiff":this.assignObjectToOutput(t,"ifd0",r);break;case"exif":this.assignObjectToOutput(t,"exif",r);break;case"xmlns":break;default:this.assignObjectToOutput(t,i,r)}else t.xmp=e}}st(Ri,"type","xmp"),st(Ri,"multiSegment",!0),Ct.set("xmp",Ri);class Zn{static findAll(t){return ss(t,/([a-zA-Z0-9-]+):([a-zA-Z0-9-]+)=("[^"]*"|'[^']*')/gm).map(Zn.unpackMatch)}static unpackMatch(t){let e=t[1],i=t[2],r=t[3].slice(1,-1);return r=as(r),new Zn(e,i,r)}constructor(t,e,i){this.ns=t,this.name=e,this.value=i}serialize(){return this.value}}class hn{static findAll(t,e,i){if(e!==void 0||i!==void 0){e=e||"[\\w\\d-]+",i=i||"[\\w\\d-]+";var r=new RegExp(`<(${e}):(${i})(#\\d+)?((\\s+?[\\w\\d-:]+=("[^"]*"|'[^']*'))*\\s*)(\\/>|>([\\s\\S]*?)<\\/\\1:\\2\\3>)`,"gm")}else r=/<([\w\d-]+):([\w\d-]+)(#\d+)?((\s+?[\w\d-:]+=("[^"]*"|'[^']*'))*\s*)(\/>|>([\s\S]*?)<\/\1:\2\3>)/gm;return ss(t,r).map(hn.unpackMatch)}static unpackMatch(t){let e=t[1],i=t[2],r=t[4],s=t[8];return new hn(e,i,r,s)}constructor(t,e,i,r){this.ns=t,this.name=e,this.attrString=i,this.innerXml=r,this.attrs=Zn.findAll(i),this.children=hn.findAll(r),this.value=this.children.length===0?as(r):void 0,this.properties=[...this.attrs,...this.children]}get isPrimitive(){return this.value!==void 0&&this.attrs.length===0&&this.children.length===0}get isListContainer(){return this.children.length===1&&this.children[0].isList}get isList(){let{ns:t,name:e}=this;return t==="rdf"&&(e==="Seq"||e==="Bag"||e==="Alt")}get isListItem(){return this.ns==="rdf"&&this.name==="li"}serialize(){if(this.properties.length===0&&this.value===void 0)return;if(this.isPrimitive)return this.value;if(this.isListContainer)return this.children[0].serialize();if(this.isList)return va(this.children.map(Sa));if(this.isListItem&&this.children.length===1&&this.attrs.length===0)return this.children[0].serialize();let t={};for(let e of this.properties)rs(e,t);return this.value!==void 0&&(t.value=this.value),$n(t)}}function rs(n,t){let e=n.serialize();e!==void 0&&(t[n.name]=e)}var Sa=n=>n.serialize(),va=n=>n.length===1?n[0]:n,Pa=(n,t)=>t[n]?t[n]:t[n]={};function ss(n,t){let e,i=[];if(!n)return i;for(;(e=t.exec(n))!==null;)i.push(e);return i}function as(n){if(function(i){return i==null||i==="null"||i==="undefined"||i===""||i.trim()===""}(n))return;let t=Number(n);if(!Number.isNaN(t))return t;let e=n.toLowerCase();return e==="true"||e!=="false"&&n.trim()}const os=["rdf:li","rdf:Seq","rdf:Bag","rdf:Alt","rdf:Description"],Ca=new RegExp(`(<|\\/)(${os.join("|")})`,"g");var ls=Object.freeze({__proto__:null,default:Ma,Exifr:De,fileParsers:xe,segmentParsers:Ct,fileReaders:_e,tagKeys:Rt,tagValues:le,tagRevivers:Ue,createDictionary:Pt,extendDictionary:an,fetchUrlAsArrayBuffer:rn,readBlobAsArrayBuffer:sn,chunkedProps:Ie,otherSegments:mn,segments:on,tiffBlocks:St,segmentsAndBlocks:Re,tiffExtractables:Ne,inheritables:gn,allFormatters:Le,Options:ln,parse:ai,gpsOnlyOptions:Qi,gps:Jr,thumbnailOnlyOptions:qi,thumbnail:Zr,thumbnailUrl:ts,orientationOnlyOptions:Ki,orientation:$i,rotations:Ji,get rotateCanvas(){return $e},get rotateCss(){return Je},rotation:es});let Fr=Kn("fs",n=>n.promises);_e.set("fs",class extends li{async readWhole(){this.chunked=!1,this.fs=await Fr;let n=await this.fs.readFile(this.input);this._swapBuffer(n)}async readChunked(){this.chunked=!0,this.fs=await Fr,await this.open(),await this.readChunk(0,this.options.firstChunkSize)}async open(){this.fh===void 0&&(this.fh=await this.fs.open(this.input,"r"),this.size=(await this.fh.stat(this.input)).size)}async _readChunk(n,t){this.fh===void 0&&await this.open(),n+t>this.size&&(t=this.size-n);var e=this.subarray(n,t,!0);return await this.fh.read(e.dataView,0,t,n),e}async close(){if(this.fh){let n=this.fh;this.fh=void 0,await n.close()}}});_e.set("base64",class extends li{constructor(...n){super(...n),this.input=this.input.replace(/^data:([^;]+);base64,/gim,""),this.size=this.input.length/4*3,this.input.endsWith("==")?this.size-=2:this.input.endsWith("=")&&(this.size-=1)}async _readChunk(n,t){let e,i,r=this.input;n===void 0?(n=0,e=0,i=0):(e=4*Math.floor(n/3),i=n-e/4*3),t===void 0&&(t=this.size);let s=n+t,o=e+4*Math.ceil(s/3);r=r.slice(e,o);let a=Math.min(t,this.size-n);if(si){let l=ri.from(r,"base64").slice(i,i+a);return this.set(l,n,!0)}{let l=this.subarray(n,a,!0),u=atob(r),d=l.toUint8();for(let c=0;c<a;c++)d[c]=u.charCodeAt(i+c);return l}}});class Ar extends oi{static canHandle(t,e){return e===18761||e===19789}extendOptions(t){let{ifd0:e,xmp:i,iptc:r,icc:s}=t;i.enabled&&e.deps.add(700),r.enabled&&e.deps.add(33723),s.enabled&&e.deps.add(34675),e.finalizeFilters()}async parse(){let{tiff:t,xmp:e,iptc:i,icc:r}=this.options;if(t.enabled||e.enabled||i.enabled||r.enabled){let s=Math.max(Di(this.options),this.options.chunkSize);await this.file.ensureChunk(0,s),this.createParser("tiff",this.file),this.parsers.tiff.parseHeader(),await this.parsers.tiff.parseIfd0Block(),this.adaptTiffPropAsSegment("xmp"),this.adaptTiffPropAsSegment("iptc"),this.adaptTiffPropAsSegment("icc")}}adaptTiffPropAsSegment(t){if(this.parsers.tiff[t]){let e=this.parsers.tiff[t];this.injectSegment(t,e)}}}st(Ar,"type","tiff"),xe.set("tiff",Ar);let ka=Kn("zlib");const Fa=["ihdr","iccp","text","itxt","exif"];class Tr extends oi{constructor(...t){super(...t),st(this,"catchError",e=>this.errors.push(e)),st(this,"metaChunks",[]),st(this,"unknownChunks",[])}static canHandle(t,e){return e===35152&&t.getUint32(0)===2303741511&&t.getUint32(4)===218765834}async parse(){let{file:t}=this;await this.findPngChunksInRange(8,t.byteLength),await this.readSegments(this.metaChunks),this.findIhdr(),this.parseTextChunks(),await this.findExif().catch(this.catchError),await this.findXmp().catch(this.catchError),await this.findIcc().catch(this.catchError)}async findPngChunksInRange(t,e){let{file:i}=this;for(;t<e;){let r=i.getUint32(t),s=i.getUint32(t+4),o=i.getString(t+4,4).toLowerCase(),a=r+4+4+4,l={type:o,offset:t,length:a,start:t+4+4,size:r,marker:s};Fa.includes(o)?this.metaChunks.push(l):this.unknownChunks.push(l),t+=a}}parseTextChunks(){let t=this.metaChunks.filter(e=>e.type==="text");for(let e of t){let[i,r]=this.file.getString(e.start,e.size).split("\0");this.injectKeyValToIhdr(i,r)}}injectKeyValToIhdr(t,e){let i=this.parsers.ihdr;i&&i.raw.set(t,e)}findIhdr(){let t=this.metaChunks.find(e=>e.type==="ihdr");t&&this.options.ihdr.enabled!==!1&&this.createParser("ihdr",t.chunk)}async findExif(){let t=this.metaChunks.find(e=>e.type==="exif");t&&this.injectSegment("tiff",t.chunk)}async findXmp(){let t=this.metaChunks.filter(e=>e.type==="itxt");for(let e of t)e.chunk.getString(0,17)==="XML:com.adobe.xmp"&&this.injectSegment("xmp",e.chunk)}async findIcc(){let t=this.metaChunks.find(a=>a.type==="iccp");if(!t)return;let{chunk:e}=t,i=e.getUint8Array(0,81),r=0;for(;r<80&&i[r]!==0;)r++;let s=r+2,o=e.getString(0,r);if(this.injectKeyValToIhdr("ProfileName",o),qn){let a=await ka,l=e.getUint8Array(s);l=a.inflateSync(l),this.injectSegment("icc",l)}}}st(Tr,"type","png"),xe.set("png",Tr),Pt(Rt,"interop",[[1,"InteropIndex"],[2,"InteropVersion"],[4096,"RelatedImageFileFormat"],[4097,"RelatedImageWidth"],[4098,"RelatedImageHeight"]]),an(Rt,"ifd0",[[11,"ProcessingSoftware"],[254,"SubfileType"],[255,"OldSubfileType"],[263,"Thresholding"],[264,"CellWidth"],[265,"CellLength"],[266,"FillOrder"],[269,"DocumentName"],[280,"MinSampleValue"],[281,"MaxSampleValue"],[285,"PageName"],[286,"XPosition"],[287,"YPosition"],[290,"GrayResponseUnit"],[297,"PageNumber"],[321,"HalftoneHints"],[322,"TileWidth"],[323,"TileLength"],[332,"InkSet"],[337,"TargetPrinter"],[18246,"Rating"],[18249,"RatingPercent"],[33550,"PixelScale"],[34264,"ModelTransform"],[34377,"PhotoshopSettings"],[50706,"DNGVersion"],[50707,"DNGBackwardVersion"],[50708,"UniqueCameraModel"],[50709,"LocalizedCameraModel"],[50736,"DNGLensInfo"],[50739,"ShadowScale"],[50740,"DNGPrivateData"],[33920,"IntergraphMatrix"],[33922,"ModelTiePoint"],[34118,"SEMInfo"],[34735,"GeoTiffDirectory"],[34736,"GeoTiffDoubleParams"],[34737,"GeoTiffAsciiParams"],[50341,"PrintIM"],[50721,"ColorMatrix1"],[50722,"ColorMatrix2"],[50723,"CameraCalibration1"],[50724,"CameraCalibration2"],[50725,"ReductionMatrix1"],[50726,"ReductionMatrix2"],[50727,"AnalogBalance"],[50728,"AsShotNeutral"],[50729,"AsShotWhiteXY"],[50730,"BaselineExposure"],[50731,"BaselineNoise"],[50732,"BaselineSharpness"],[50734,"LinearResponseLimit"],[50735,"CameraSerialNumber"],[50741,"MakerNoteSafety"],[50778,"CalibrationIlluminant1"],[50779,"CalibrationIlluminant2"],[50781,"RawDataUniqueID"],[50827,"OriginalRawFileName"],[50828,"OriginalRawFileData"],[50831,"AsShotICCProfile"],[50832,"AsShotPreProfileMatrix"],[50833,"CurrentICCProfile"],[50834,"CurrentPreProfileMatrix"],[50879,"ColorimetricReference"],[50885,"SRawType"],[50898,"PanasonicTitle"],[50899,"PanasonicTitle2"],[50931,"CameraCalibrationSig"],[50932,"ProfileCalibrationSig"],[50933,"ProfileIFD"],[50934,"AsShotProfileName"],[50936,"ProfileName"],[50937,"ProfileHueSatMapDims"],[50938,"ProfileHueSatMapData1"],[50939,"ProfileHueSatMapData2"],[50940,"ProfileToneCurve"],[50941,"ProfileEmbedPolicy"],[50942,"ProfileCopyright"],[50964,"ForwardMatrix1"],[50965,"ForwardMatrix2"],[50966,"PreviewApplicationName"],[50967,"PreviewApplicationVersion"],[50968,"PreviewSettingsName"],[50969,"PreviewSettingsDigest"],[50970,"PreviewColorSpace"],[50971,"PreviewDateTime"],[50972,"RawImageDigest"],[50973,"OriginalRawFileDigest"],[50981,"ProfileLookTableDims"],[50982,"ProfileLookTableData"],[51043,"TimeCodes"],[51044,"FrameRate"],[51058,"TStop"],[51081,"ReelName"],[51089,"OriginalDefaultFinalSize"],[51090,"OriginalBestQualitySize"],[51091,"OriginalDefaultCropSize"],[51105,"CameraLabel"],[51107,"ProfileHueSatMapEncoding"],[51108,"ProfileLookTableEncoding"],[51109,"BaselineExposureOffset"],[51110,"DefaultBlackRender"],[51111,"NewRawImageDigest"],[51112,"RawToPreviewGain"]]);let Ir=[[273,"StripOffsets"],[279,"StripByteCounts"],[288,"FreeOffsets"],[289,"FreeByteCounts"],[291,"GrayResponseCurve"],[292,"T4Options"],[293,"T6Options"],[300,"ColorResponseUnit"],[320,"ColorMap"],[324,"TileOffsets"],[325,"TileByteCounts"],[326,"BadFaxLines"],[327,"CleanFaxData"],[328,"ConsecutiveBadFaxLines"],[330,"SubIFD"],[333,"InkNames"],[334,"NumberofInks"],[336,"DotRange"],[338,"ExtraSamples"],[339,"SampleFormat"],[340,"SMinSampleValue"],[341,"SMaxSampleValue"],[342,"TransferRange"],[343,"ClipPath"],[344,"XClipPathUnits"],[345,"YClipPathUnits"],[346,"Indexed"],[347,"JPEGTables"],[351,"OPIProxy"],[400,"GlobalParametersIFD"],[401,"ProfileType"],[402,"FaxProfile"],[403,"CodingMethods"],[404,"VersionYear"],[405,"ModeNumber"],[433,"Decode"],[434,"DefaultImageColor"],[435,"T82Options"],[437,"JPEGTables"],[512,"JPEGProc"],[515,"JPEGRestartInterval"],[517,"JPEGLosslessPredictors"],[518,"JPEGPointTransforms"],[519,"JPEGQTables"],[520,"JPEGDCTables"],[521,"JPEGACTables"],[559,"StripRowCounts"],[999,"USPTOMiscellaneous"],[18247,"XP_DIP_XML"],[18248,"StitchInfo"],[28672,"SonyRawFileType"],[28688,"SonyToneCurve"],[28721,"VignettingCorrection"],[28722,"VignettingCorrParams"],[28724,"ChromaticAberrationCorrection"],[28725,"ChromaticAberrationCorrParams"],[28726,"DistortionCorrection"],[28727,"DistortionCorrParams"],[29895,"SonyCropTopLeft"],[29896,"SonyCropSize"],[32781,"ImageID"],[32931,"WangTag1"],[32932,"WangAnnotation"],[32933,"WangTag3"],[32934,"WangTag4"],[32953,"ImageReferencePoints"],[32954,"RegionXformTackPoint"],[32955,"WarpQuadrilateral"],[32956,"AffineTransformMat"],[32995,"Matteing"],[32996,"DataType"],[32997,"ImageDepth"],[32998,"TileDepth"],[33300,"ImageFullWidth"],[33301,"ImageFullHeight"],[33302,"TextureFormat"],[33303,"WrapModes"],[33304,"FovCot"],[33305,"MatrixWorldToScreen"],[33306,"MatrixWorldToCamera"],[33405,"Model2"],[33421,"CFARepeatPatternDim"],[33422,"CFAPattern2"],[33423,"BatteryLevel"],[33424,"KodakIFD"],[33445,"MDFileTag"],[33446,"MDScalePixel"],[33447,"MDColorTable"],[33448,"MDLabName"],[33449,"MDSampleInfo"],[33450,"MDPrepDate"],[33451,"MDPrepTime"],[33452,"MDFileUnits"],[33589,"AdventScale"],[33590,"AdventRevision"],[33628,"UIC1Tag"],[33629,"UIC2Tag"],[33630,"UIC3Tag"],[33631,"UIC4Tag"],[33918,"IntergraphPacketData"],[33919,"IntergraphFlagRegisters"],[33921,"INGRReserved"],[34016,"Site"],[34017,"ColorSequence"],[34018,"IT8Header"],[34019,"RasterPadding"],[34020,"BitsPerRunLength"],[34021,"BitsPerExtendedRunLength"],[34022,"ColorTable"],[34023,"ImageColorIndicator"],[34024,"BackgroundColorIndicator"],[34025,"ImageColorValue"],[34026,"BackgroundColorValue"],[34027,"PixelIntensityRange"],[34028,"TransparencyIndicator"],[34029,"ColorCharacterization"],[34030,"HCUsage"],[34031,"TrapIndicator"],[34032,"CMYKEquivalent"],[34152,"AFCP_IPTC"],[34232,"PixelMagicJBIGOptions"],[34263,"JPLCartoIFD"],[34306,"WB_GRGBLevels"],[34310,"LeafData"],[34687,"TIFF_FXExtensions"],[34688,"MultiProfiles"],[34689,"SharedData"],[34690,"T88Options"],[34732,"ImageLayer"],[34750,"JBIGOptions"],[34856,"Opto-ElectricConvFactor"],[34857,"Interlace"],[34908,"FaxRecvParams"],[34909,"FaxSubAddress"],[34910,"FaxRecvTime"],[34929,"FedexEDR"],[34954,"LeafSubIFD"],[37387,"FlashEnergy"],[37388,"SpatialFrequencyResponse"],[37389,"Noise"],[37390,"FocalPlaneXResolution"],[37391,"FocalPlaneYResolution"],[37392,"FocalPlaneResolutionUnit"],[37397,"ExposureIndex"],[37398,"TIFF-EPStandardID"],[37399,"SensingMethod"],[37434,"CIP3DataFile"],[37435,"CIP3Sheet"],[37436,"CIP3Side"],[37439,"StoNits"],[37679,"MSDocumentText"],[37680,"MSPropertySetStorage"],[37681,"MSDocumentTextPosition"],[37724,"ImageSourceData"],[40965,"InteropIFD"],[40976,"SamsungRawPointersOffset"],[40977,"SamsungRawPointersLength"],[41217,"SamsungRawByteOrder"],[41218,"SamsungRawUnknown"],[41484,"SpatialFrequencyResponse"],[41485,"Noise"],[41489,"ImageNumber"],[41490,"SecurityClassification"],[41491,"ImageHistory"],[41494,"TIFF-EPStandardID"],[41995,"DeviceSettingDescription"],[42112,"GDALMetadata"],[42113,"GDALNoData"],[44992,"ExpandSoftware"],[44993,"ExpandLens"],[44994,"ExpandFilm"],[44995,"ExpandFilterLens"],[44996,"ExpandScanner"],[44997,"ExpandFlashLamp"],[46275,"HasselbladRawImage"],[48129,"PixelFormat"],[48130,"Transformation"],[48131,"Uncompressed"],[48132,"ImageType"],[48256,"ImageWidth"],[48257,"ImageHeight"],[48258,"WidthResolution"],[48259,"HeightResolution"],[48320,"ImageOffset"],[48321,"ImageByteCount"],[48322,"AlphaOffset"],[48323,"AlphaByteCount"],[48324,"ImageDataDiscard"],[48325,"AlphaDataDiscard"],[50215,"OceScanjobDesc"],[50216,"OceApplicationSelector"],[50217,"OceIDNumber"],[50218,"OceImageLogic"],[50255,"Annotations"],[50459,"HasselbladExif"],[50547,"OriginalFileName"],[50560,"USPTOOriginalContentType"],[50656,"CR2CFAPattern"],[50710,"CFAPlaneColor"],[50711,"CFALayout"],[50712,"LinearizationTable"],[50713,"BlackLevelRepeatDim"],[50714,"BlackLevel"],[50715,"BlackLevelDeltaH"],[50716,"BlackLevelDeltaV"],[50717,"WhiteLevel"],[50718,"DefaultScale"],[50719,"DefaultCropOrigin"],[50720,"DefaultCropSize"],[50733,"BayerGreenSplit"],[50737,"ChromaBlurRadius"],[50738,"AntiAliasStrength"],[50752,"RawImageSegmentation"],[50780,"BestQualityScale"],[50784,"AliasLayerMetadata"],[50829,"ActiveArea"],[50830,"MaskedAreas"],[50935,"NoiseReductionApplied"],[50974,"SubTileBlockSize"],[50975,"RowInterleaveFactor"],[51008,"OpcodeList1"],[51009,"OpcodeList2"],[51022,"OpcodeList3"],[51041,"NoiseProfile"],[51114,"CacheVersion"],[51125,"DefaultUserCrop"],[51157,"NikonNEFInfo"],[65024,"KdcIFD"]];an(Rt,"ifd0",Ir),an(Rt,"exif",Ir),Pt(le,"gps",[[23,{M:"Magnetic North",T:"True North"}],[25,{K:"Kilometers",M:"Miles",N:"Nautical Miles"}]]);class Ni extends oe{static canHandle(t,e){return t.getUint8(e+1)===224&&t.getUint32(e+4)===1246120262&&t.getUint8(e+8)===0}parse(){return this.parseTags(),this.translate(),this.output}parseTags(){this.raw=new Map([[0,this.chunk.getUint16(0)],[2,this.chunk.getUint8(2)],[3,this.chunk.getUint16(3)],[5,this.chunk.getUint16(5)],[7,this.chunk.getUint8(7)],[8,this.chunk.getUint8(8)]])}}st(Ni,"type","jfif"),st(Ni,"headerLength",9),Ct.set("jfif",Ni),Pt(Rt,"jfif",[[0,"JFIFVersion"],[2,"ResolutionUnit"],[3,"XResolution"],[5,"YResolution"],[7,"ThumbnailWidth"],[8,"ThumbnailHeight"]]);class Rr extends oe{parse(){return this.parseTags(),this.translate(),this.output}parseTags(){this.raw=new Map([[0,this.chunk.getUint32(0)],[4,this.chunk.getUint32(4)],[8,this.chunk.getUint8(8)],[9,this.chunk.getUint8(9)],[10,this.chunk.getUint8(10)],[11,this.chunk.getUint8(11)],[12,this.chunk.getUint8(12)],...Array.from(this.raw)])}}st(Rr,"type","ihdr"),Ct.set("ihdr",Rr),Pt(Rt,"ihdr",[[0,"ImageWidth"],[4,"ImageHeight"],[8,"BitDepth"],[9,"ColorType"],[10,"Compression"],[11,"Filter"],[12,"Interlace"]]),Pt(le,"ihdr",[[9,{0:"Grayscale",2:"RGB",3:"Palette",4:"Grayscale with Alpha",6:"RGB with Alpha",DEFAULT:"Unknown"}],[10,{0:"Deflate/Inflate",DEFAULT:"Unknown"}],[11,{0:"Adaptive",DEFAULT:"Unknown"}],[12,{0:"Noninterlaced",1:"Adam7 Interlace",DEFAULT:"Unknown"}]]);class Qn extends oe{static canHandle(t,e){return t.getUint8(e+1)===226&&t.getUint32(e+4)===1229144927}static findPosition(t,e){let i=super.findPosition(t,e);return i.chunkNumber=t.getUint8(e+16),i.chunkCount=t.getUint8(e+17),i.multiSegment=i.chunkCount>1,i}static handleMultiSegments(t){return function(e){let i=function(r){let s=r[0].constructor,o=0;for(let u of r)o+=u.length;let a=new s(o),l=0;for(let u of r)a.set(u,l),l+=u.length;return a}(e.map(r=>r.chunk.toUint8()));return new jt(i)}(t)}parse(){return this.raw=new Map,this.parseHeader(),this.parseTags(),this.translate(),this.output}parseHeader(){let{raw:t}=this;this.chunk.byteLength<84&&kt("ICC header is too short");for(let[e,i]of Object.entries(Aa)){e=parseInt(e,10);let r=i(this.chunk,e);r!=="\0\0\0\0"&&t.set(e,r)}}parseTags(){let t,e,i,r,s,{raw:o}=this,a=this.chunk.getUint32(128),l=132,u=this.chunk.byteLength;for(;a--;){if(t=this.chunk.getString(l,4),e=this.chunk.getUint32(l+4),i=this.chunk.getUint32(l+8),r=this.chunk.getString(e,4),e+i>u)return void console.warn("reached the end of the first ICC chunk. Enable options.tiff.multiSegment to read all ICC segments.");s=this.parseTag(r,e,i),s!==void 0&&s!=="\0\0\0\0"&&o.set(t,s),l+=12}}parseTag(t,e,i){switch(t){case"desc":return this.parseDesc(e);case"mluc":return this.parseMluc(e);case"text":return this.parseText(e,i);case"sig ":return this.parseSig(e)}if(!(e+i>this.chunk.byteLength))return this.chunk.getUint8Array(e,i)}parseDesc(t){let e=this.chunk.getUint32(t+8)-1;return Ke(this.chunk.getString(t+12,e))}parseText(t,e){return Ke(this.chunk.getString(t+8,e-8))}parseSig(t){return Ke(this.chunk.getString(t+8,4))}parseMluc(t){let{chunk:e}=this,i=e.getUint32(t+8),r=e.getUint32(t+12),s=t+16,o=[];for(let a=0;a<i;a++){let l=e.getString(s+0,2),u=e.getString(s+2,2),d=e.getUint32(s+4),c=e.getUint32(s+8)+t,m=Ke(e.getUnicodeString(c,d));o.push({lang:l,country:u,text:m}),s+=r}return i===1?o[0].text:o}translateValue(t,e){return typeof t=="string"?e[t]||e[t.toLowerCase()]||t:e[t]||t}}st(Qn,"type","icc"),st(Qn,"multiSegment",!0),st(Qn,"headerLength",18);const Aa={4:Se,8:function(n,t){return[n.getUint8(t),n.getUint8(t+1)>>4,n.getUint8(t+1)%16].map(e=>e.toString(10)).join(".")},12:Se,16:Se,20:Se,24:function(n,t){const e=n.getUint16(t),i=n.getUint16(t+2)-1,r=n.getUint16(t+4),s=n.getUint16(t+6),o=n.getUint16(t+8),a=n.getUint16(t+10);return new Date(Date.UTC(e,i,r,s,o,a))},36:Se,40:Se,48:Se,52:Se,64:(n,t)=>n.getUint32(t),80:Se};function Se(n,t){return Ke(n.getString(t,4))}Ct.set("icc",Qn),Pt(Rt,"icc",[[4,"ProfileCMMType"],[8,"ProfileVersion"],[12,"ProfileClass"],[16,"ColorSpaceData"],[20,"ProfileConnectionSpace"],[24,"ProfileDateTime"],[36,"ProfileFileSignature"],[40,"PrimaryPlatform"],[44,"CMMFlags"],[48,"DeviceManufacturer"],[52,"DeviceModel"],[56,"DeviceAttributes"],[64,"RenderingIntent"],[68,"ConnectionSpaceIlluminant"],[80,"ProfileCreator"],[84,"ProfileID"],["Header","ProfileHeader"],["MS00","WCSProfiles"],["bTRC","BlueTRC"],["bXYZ","BlueMatrixColumn"],["bfd","UCRBG"],["bkpt","MediaBlackPoint"],["calt","CalibrationDateTime"],["chad","ChromaticAdaptation"],["chrm","Chromaticity"],["ciis","ColorimetricIntentImageState"],["clot","ColorantTableOut"],["clro","ColorantOrder"],["clrt","ColorantTable"],["cprt","ProfileCopyright"],["crdi","CRDInfo"],["desc","ProfileDescription"],["devs","DeviceSettings"],["dmdd","DeviceModelDesc"],["dmnd","DeviceMfgDesc"],["dscm","ProfileDescriptionML"],["fpce","FocalPlaneColorimetryEstimates"],["gTRC","GreenTRC"],["gXYZ","GreenMatrixColumn"],["gamt","Gamut"],["kTRC","GrayTRC"],["lumi","Luminance"],["meas","Measurement"],["meta","Metadata"],["mmod","MakeAndModel"],["ncl2","NamedColor2"],["ncol","NamedColor"],["ndin","NativeDisplayInfo"],["pre0","Preview0"],["pre1","Preview1"],["pre2","Preview2"],["ps2i","PS2RenderingIntent"],["ps2s","PostScript2CSA"],["psd0","PostScript2CRD0"],["psd1","PostScript2CRD1"],["psd2","PostScript2CRD2"],["psd3","PostScript2CRD3"],["pseq","ProfileSequenceDesc"],["psid","ProfileSequenceIdentifier"],["psvm","PS2CRDVMSize"],["rTRC","RedTRC"],["rXYZ","RedMatrixColumn"],["resp","OutputResponse"],["rhoc","ReflectionHardcopyOrigColorimetry"],["rig0","PerceptualRenderingIntentGamut"],["rig2","SaturationRenderingIntentGamut"],["rpoc","ReflectionPrintOutputColorimetry"],["sape","SceneAppearanceEstimates"],["scoe","SceneColorimetryEstimates"],["scrd","ScreeningDesc"],["scrn","Screening"],["targ","CharTarget"],["tech","Technology"],["vcgt","VideoCardGamma"],["view","ViewingConditions"],["vued","ViewingCondDesc"],["wtpt","MediaWhitePoint"]]);const Gn={"4d2p":"Erdt Systems",AAMA:"Aamazing Technologies",ACER:"Acer",ACLT:"Acolyte Color Research",ACTI:"Actix Sytems",ADAR:"Adara Technology",ADBE:"Adobe",ADI:"ADI Systems",AGFA:"Agfa Graphics",ALMD:"Alps Electric",ALPS:"Alps Electric",ALWN:"Alwan Color Expertise",AMTI:"Amiable Technologies",AOC:"AOC International",APAG:"Apago",APPL:"Apple Computer",AST:"AST","AT&T":"AT&T",BAEL:"BARBIERI electronic",BRCO:"Barco NV",BRKP:"Breakpoint",BROT:"Brother",BULL:"Bull",BUS:"Bus Computer Systems","C-IT":"C-Itoh",CAMR:"Intel",CANO:"Canon",CARR:"Carroll Touch",CASI:"Casio",CBUS:"Colorbus PL",CEL:"Crossfield",CELx:"Crossfield",CGS:"CGS Publishing Technologies International",CHM:"Rochester Robotics",CIGL:"Colour Imaging Group, London",CITI:"Citizen",CL00:"Candela",CLIQ:"Color IQ",CMCO:"Chromaco",CMiX:"CHROMiX",COLO:"Colorgraphic Communications",COMP:"Compaq",COMp:"Compeq/Focus Technology",CONR:"Conrac Display Products",CORD:"Cordata Technologies",CPQ:"Compaq",CPRO:"ColorPro",CRN:"Cornerstone",CTX:"CTX International",CVIS:"ColorVision",CWC:"Fujitsu Laboratories",DARI:"Darius Technology",DATA:"Dataproducts",DCP:"Dry Creek Photo",DCRC:"Digital Contents Resource Center, Chung-Ang University",DELL:"Dell Computer",DIC:"Dainippon Ink and Chemicals",DICO:"Diconix",DIGI:"Digital","DL&C":"Digital Light & Color",DPLG:"Doppelganger",DS:"Dainippon Screen",DSOL:"DOOSOL",DUPN:"DuPont",EPSO:"Epson",ESKO:"Esko-Graphics",ETRI:"Electronics and Telecommunications Research Institute",EVER:"Everex Systems",EXAC:"ExactCODE",Eizo:"Eizo",FALC:"Falco Data Products",FF:"Fuji Photo Film",FFEI:"FujiFilm Electronic Imaging",FNRD:"Fnord Software",FORA:"Fora",FORE:"Forefront Technology",FP:"Fujitsu",FPA:"WayTech Development",FUJI:"Fujitsu",FX:"Fuji Xerox",GCC:"GCC Technologies",GGSL:"Global Graphics Software",GMB:"Gretagmacbeth",GMG:"GMG",GOLD:"GoldStar Technology",GOOG:"Google",GPRT:"Giantprint",GTMB:"Gretagmacbeth",GVC:"WayTech Development",GW2K:"Sony",HCI:"HCI",HDM:"Heidelberger Druckmaschinen",HERM:"Hermes",HITA:"Hitachi America",HP:"Hewlett-Packard",HTC:"Hitachi",HiTi:"HiTi Digital",IBM:"IBM",IDNT:"Scitex",IEC:"Hewlett-Packard",IIYA:"Iiyama North America",IKEG:"Ikegami Electronics",IMAG:"Image Systems",IMI:"Ingram Micro",INTC:"Intel",INTL:"N/A (INTL)",INTR:"Intra Electronics",IOCO:"Iocomm International Technology",IPS:"InfoPrint Solutions Company",IRIS:"Scitex",ISL:"Ichikawa Soft Laboratory",ITNL:"N/A (ITNL)",IVM:"IVM",IWAT:"Iwatsu Electric",Idnt:"Scitex",Inca:"Inca Digital Printers",Iris:"Scitex",JPEG:"Joint Photographic Experts Group",JSFT:"Jetsoft Development",JVC:"JVC Information Products",KART:"Scitex",KFC:"KFC Computek Components",KLH:"KLH Computers",KMHD:"Konica Minolta",KNCA:"Konica",KODA:"Kodak",KYOC:"Kyocera",Kart:"Scitex",LCAG:"Leica",LCCD:"Leeds Colour",LDAK:"Left Dakota",LEAD:"Leading Technology",LEXM:"Lexmark International",LINK:"Link Computer",LINO:"Linotronic",LITE:"Lite-On",Leaf:"Leaf",Lino:"Linotronic",MAGC:"Mag Computronic",MAGI:"MAG Innovision",MANN:"Mannesmann",MICN:"Micron Technology",MICR:"Microtek",MICV:"Microvitec",MINO:"Minolta",MITS:"Mitsubishi Electronics America",MITs:"Mitsuba",MNLT:"Minolta",MODG:"Modgraph",MONI:"Monitronix",MONS:"Monaco Systems",MORS:"Morse Technology",MOTI:"Motive Systems",MSFT:"Microsoft",MUTO:"MUTOH INDUSTRIES",Mits:"Mitsubishi Electric",NANA:"NANAO",NEC:"NEC",NEXP:"NexPress Solutions",NISS:"Nissei Sangyo America",NKON:"Nikon",NONE:"none",OCE:"Oce Technologies",OCEC:"OceColor",OKI:"Oki",OKID:"Okidata",OKIP:"Okidata",OLIV:"Olivetti",OLYM:"Olympus",ONYX:"Onyx Graphics",OPTI:"Optiquest",PACK:"Packard Bell",PANA:"Matsushita Electric Industrial",PANT:"Pantone",PBN:"Packard Bell",PFU:"PFU",PHIL:"Philips Consumer Electronics",PNTX:"HOYA",POne:"Phase One A/S",PREM:"Premier Computer Innovations",PRIN:"Princeton Graphic Systems",PRIP:"Princeton Publishing Labs",QLUX:"Hong Kong",QMS:"QMS",QPCD:"QPcard AB",QUAD:"QuadLaser",QUME:"Qume",RADI:"Radius",RDDx:"Integrated Color Solutions",RDG:"Roland DG",REDM:"REDMS Group",RELI:"Relisys",RGMS:"Rolf Gierling Multitools",RICO:"Ricoh",RNLD:"Edmund Ronald",ROYA:"Royal",RPC:"Ricoh Printing Systems",RTL:"Royal Information Electronics",SAMP:"Sampo",SAMS:"Samsung",SANT:"Jaime Santana Pomares",SCIT:"Scitex",SCRN:"Dainippon Screen",SDP:"Scitex",SEC:"Samsung",SEIK:"Seiko Instruments",SEIk:"Seikosha",SGUY:"ScanGuy.com",SHAR:"Sharp Laboratories",SICC:"International Color Consortium",SONY:"Sony",SPCL:"SpectraCal",STAR:"Star",STC:"Sampo Technology",Scit:"Scitex",Sdp:"Scitex",Sony:"Sony",TALO:"Talon Technology",TAND:"Tandy",TATU:"Tatung",TAXA:"TAXAN America",TDS:"Tokyo Denshi Sekei",TECO:"TECO Information Systems",TEGR:"Tegra",TEKT:"Tektronix",TI:"Texas Instruments",TMKR:"TypeMaker",TOSB:"Toshiba",TOSH:"Toshiba",TOTK:"TOTOKU ELECTRIC",TRIU:"Triumph",TSBT:"Toshiba",TTX:"TTX Computer Products",TVM:"TVM Professional Monitor",TW:"TW Casper",ULSX:"Ulead Systems",UNIS:"Unisys",UTZF:"Utz Fehlau & Sohn",VARI:"Varityper",VIEW:"Viewsonic",VISL:"Visual communication",VIVO:"Vivo Mobile Communication",WANG:"Wang",WLBR:"Wilbur Imaging",WTG2:"Ware To Go",WYSE:"WYSE Technology",XERX:"Xerox",XRIT:"X-Rite",ZRAN:"Zoran",Zebr:"Zebra Technologies",appl:"Apple Computer",bICC:"basICColor",berg:"bergdesign",ceyd:"Integrated Color Solutions",clsp:"MacDermid ColorSpan",ds:"Dainippon Screen",dupn:"DuPont",ffei:"FujiFilm Electronic Imaging",flux:"FluxData",iris:"Scitex",kart:"Scitex",lcms:"Little CMS",lino:"Linotronic",none:"none",ob4d:"Erdt Systems",obic:"Medigraph",quby:"Qubyx Sarl",scit:"Scitex",scrn:"Dainippon Screen",sdp:"Scitex",siwi:"SIWI GRAFIKA",yxym:"YxyMaster"},Nr={scnr:"Scanner",mntr:"Monitor",prtr:"Printer",link:"Device Link",abst:"Abstract",spac:"Color Space Conversion Profile",nmcl:"Named Color",cenc:"ColorEncodingSpace profile",mid:"MultiplexIdentification profile",mlnk:"MultiplexLink profile",mvis:"MultiplexVisualization profile",nkpf:"Nikon Input Device Profile (NON-STANDARD!)"};Pt(le,"icc",[[4,Gn],[12,Nr],[40,Object.assign({},Gn,Nr)],[48,Gn],[80,Gn],[64,{0:"Perceptual",1:"Relative Colorimetric",2:"Saturation",3:"Absolute Colorimetric"}],["tech",{amd:"Active Matrix Display",crt:"Cathode Ray Tube Display",kpcd:"Photo CD",pmd:"Passive Matrix Display",dcam:"Digital Camera",dcpj:"Digital Cinema Projector",dmpc:"Digital Motion Picture Camera",dsub:"Dye Sublimation Printer",epho:"Electrophotographic Printer",esta:"Electrostatic Printer",flex:"Flexography",fprn:"Film Writer",fscn:"Film Scanner",grav:"Gravure",ijet:"Ink Jet Printer",imgs:"Photo Image Setter",mpfr:"Motion Picture Film Recorder",mpfs:"Motion Picture Film Scanner",offs:"Offset Lithography",pjtv:"Projection Television",rpho:"Photographic Paper Printer",rscn:"Reflective Scanner",silk:"Silkscreen",twax:"Thermal Wax Printer",vidc:"Video Camera",vidm:"Video Monitor"}]]);class Xn extends oe{static canHandle(t,e,i){return t.getUint8(e+1)===237&&t.getString(e+4,9)==="Photoshop"&&this.containsIptc8bim(t,e,i)!==void 0}static headerLength(t,e,i){let r,s=this.containsIptc8bim(t,e,i);if(s!==void 0)return r=t.getUint8(e+s+7),r%2!=0&&(r+=1),r===0&&(r=4),s+8+r}static containsIptc8bim(t,e,i){for(let r=0;r<i;r++)if(this.isIptcSegmentHead(t,e+r))return r}static isIptcSegmentHead(t,e){return t.getUint8(e)===56&&t.getUint32(e)===943868237&&t.getUint16(e+4)===1028}parse(){let{raw:t}=this,e=this.chunk.byteLength-1,i=!1;for(let r=0;r<e;r++)if(this.chunk.getUint8(r)===28&&this.chunk.getUint8(r+1)===2){i=!0;let s=this.chunk.getUint16(r+3),o=this.chunk.getUint8(r+2),a=this.chunk.getLatin1String(r+5,s);t.set(o,this.pluralizeValue(t.get(o),a)),r+=4+s}else if(i)break;return this.translate(),this.output}pluralizeValue(t,e){return t!==void 0?t instanceof Array?(t.push(e),t):[t,e]:e}}st(Xn,"type","iptc"),st(Xn,"translateValues",!1),st(Xn,"reviveValues",!1),Ct.set("iptc",Xn),Pt(Rt,"iptc",[[0,"ApplicationRecordVersion"],[3,"ObjectTypeReference"],[4,"ObjectAttributeReference"],[5,"ObjectName"],[7,"EditStatus"],[8,"EditorialUpdate"],[10,"Urgency"],[12,"SubjectReference"],[15,"Category"],[20,"SupplementalCategories"],[22,"FixtureIdentifier"],[25,"Keywords"],[26,"ContentLocationCode"],[27,"ContentLocationName"],[30,"ReleaseDate"],[35,"ReleaseTime"],[37,"ExpirationDate"],[38,"ExpirationTime"],[40,"SpecialInstructions"],[42,"ActionAdvised"],[45,"ReferenceService"],[47,"ReferenceDate"],[50,"ReferenceNumber"],[55,"DateCreated"],[60,"TimeCreated"],[62,"DigitalCreationDate"],[63,"DigitalCreationTime"],[65,"OriginatingProgram"],[70,"ProgramVersion"],[75,"ObjectCycle"],[80,"Byline"],[85,"BylineTitle"],[90,"City"],[92,"Sublocation"],[95,"State"],[100,"CountryCode"],[101,"Country"],[103,"OriginalTransmissionReference"],[105,"Headline"],[110,"Credit"],[115,"Source"],[116,"CopyrightNotice"],[118,"Contact"],[120,"Caption"],[121,"LocalCaption"],[122,"Writer"],[125,"RasterizedCaption"],[130,"ImageType"],[131,"ImageOrientation"],[135,"LanguageIdentifier"],[150,"AudioType"],[151,"AudioSamplingRate"],[152,"AudioSamplingResolution"],[153,"AudioDuration"],[154,"AudioOutcue"],[184,"JobID"],[185,"MasterDocumentID"],[186,"ShortDocumentID"],[187,"UniqueDocumentID"],[188,"OwnerID"],[200,"ObjectPreviewFileFormat"],[201,"ObjectPreviewFileVersion"],[202,"ObjectPreviewData"],[221,"Prefs"],[225,"ClassifyState"],[228,"SimilarityIndex"],[230,"DocumentNotes"],[231,"DocumentHistory"],[232,"ExifCameraInfo"],[255,"CatalogSets"]]),Pt(le,"iptc",[[10,{0:"0 (reserved)",1:"1 (most urgent)",2:"2",3:"3",4:"4",5:"5 (normal urgency)",6:"6",7:"7",8:"8 (least urgent)",9:"9 (user-defined priority)"}],[75,{a:"Morning",b:"Both Morning and Evening",p:"Evening"}],[131,{L:"Landscape",P:"Portrait",S:"Square"}]]);let Li=null;async function cs(){if(!Li)try{const i=(await import("./joraw-1Lq5hXK7.js")).default;if(typeof i!="function")throw new Error("JoRaw WASM import failed");const r=new URL("/assets/joraw-DraTMNgX.wasm",import.meta.url).href;Li=i({locateFile:(s,o)=>s.endsWith("joraw.wasm")?r:o+s})}catch(e){throw console.error("Failed to load joraw.js:",e),e}const n=await Li,t=n.LibRaw||n.JoRaw;if(!t)throw new Error("JoRaw class not found");return t}function Ta(n,t,e){const i=e.type===1||e.type===2||e.type===7?1:e.type===3?2:e.type===4||e.type===9||e.type===11||e.type===13?4:e.type===5||e.type===10||e.type===12?8:0;if(i===0)throw new Error(`Unsupported TIFF field type: ${e.type}`);const s=e.count*i<=4?e.valueFieldOffset:n.getUint32(e.valueFieldOffset,t),o=[];for(let a=0;a<e.count;a++){const l=s+a*i;if(l<0||l+i>n.byteLength)throw new Error("Invalid TIFF field offset");if(e.type===1)o.push(n.getUint8(l));else if(e.type===2||e.type===7)o.push(n.getUint8(l));else if(e.type===3)o.push(n.getUint16(l,t));else if(e.type===5){const u=n.getUint32(l+4,t);o.push(u?n.getUint32(l,t)/u:0)}else if(e.type===13)o.push(n.getUint32(l,t));else if(e.type===9)o.push(n.getInt32(l,t));else if(e.type===10){const u=n.getInt32(l+4,t);o.push(u?n.getInt32(l,t)/u:0)}else e.type===11?o.push(n.getFloat32(l,t)):e.type===12?o.push(n.getFloat64(l,t)):o.push(n.getUint32(l,t))}return o}function Lr(n,t,e){if(!Number.isFinite(e)||e<=0||e+2>n.byteLength)throw new Error("Invalid TIFF IFD offset");const i=n.getUint16(e,t);if(e+2+i*12+4>n.byteLength)throw new Error("Corrupt TIFF IFD");const s=new Map;for(let o=0;o<i;o++){const a=e+2+o*12,l=n.getUint16(a,t),u=n.getUint16(a+2,t),d=n.getUint32(a+4,t);s.set(l,{type:u,count:d,valueFieldOffset:a+8})}return s}function Wt(n,t,e,i){const r=e.get(i);return r?Ta(n,t,r):[]}const us=(...n)=>{for(const t of n){if(Array.isArray(t)||ArrayBuffer.isView(t)){const i=Array.from(t).map(Number).filter(r=>Number.isFinite(r)&&r>0);if(i.length)return Math.max(...i);continue}const e=Number(t);if(Number.isFinite(e)&&e>0)return e}return null},Ia=n=>{const e=(Array.isArray(n)||ArrayBuffer.isView(n)?Array.from(n).map(Number):[Number(n)]).map(i=>Number.isFinite(i)?i:0);return e.length>=3?[e[0]||0,e[1]||0,e[1]||0,e[2]||0]:e.length===1?[e[0]||0,e[0]||0,e[0]||0,e[0]||0]:[0,0,0,0]};function Ra(n,t={}){var T,C,k,M,S,R;if(n.byteLength<8)return null;const e=new DataView(n.buffer,n.byteOffset,n.byteLength),i=e.getUint16(0,!1),r=i===18761;if(!r&&i!==19789||e.getUint16(2,r)!==42)return null;const s=e.getUint32(4,r);let o;try{o=Lr(e,r,s)}catch{return null}const a=new Set([s]);for(const A of Wt(e,r,o,330))A>0&&a.add(A);let l=null;for(const A of a)try{const L=A===s?o:Lr(e,r,A),X=Wt(e,r,L,256)[0]||0,D=Wt(e,r,L,257)[0]||0,E=Wt(e,r,L,258),U=Wt(e,r,L,277)[0]||E.length||1,O=E.length===1?new Array(U).fill(E[0]):E,z=Wt(e,r,L,259)[0]||1,J=Wt(e,r,L,262)[0]||0,V=Wt(e,r,L,284)[0]||1,Q=Wt(e,r,L,339),W=Q.length===0?new Array(U).fill(1):Q.length===1?new Array(U).fill(Q[0]):Q,tt=Wt(e,r,L,273),at=Wt(e,r,L,279),B=Wt(e,r,L,278)[0]||D,Z=W.slice(0,U).every(I=>I===0||I===1);if(!(J===34892&&z===1&&V===1&&X>0&&D>0&&U>=3&&O.length>=U&&O.slice(0,U).every(I=>I===16)&&Z&&tt.length>0&&tt.length===at.length))continue;(!l||X*D>l.width*l.height)&&(l={entries:L,width:X,height:D,samplesPerPixel:U,rowsPerStrip:B,stripOffsets:tt,stripByteCounts:at,bitsPerSample:O,photometric:J})}catch{continue}if(!l)return null;const{entries:u,width:d,height:c,samplesPerPixel:m,rowsPerStrip:h,stripOffsets:f,stripByteCounts:p,bitsPerSample:g,photometric:x}=l,y=d*c;if(!Number.isSafeInteger(y)||y<=0)return null;const b=new Uint16Array(y*3);for(let A=0;A<f.length;A++){const L=A*h;if(L>=c)break;const X=Math.min(h,c-L),D=X*d*m*2,E=f[A],U=p[A];if(E<0||U<D||E+D>n.byteLength)throw new Error("Invalid LinearRaw DNG strip bounds");const O=L*d*3;if(m===3&&r&&!(n.byteOffset+E&1)){b.set(new Uint16Array(n.buffer,n.byteOffset+E,X*d*3),O);continue}const z=new DataView(n.buffer,n.byteOffset+E,D);let J=0,V=O;for(let Q=0;Q<X*d;Q++)b[V++]=z.getUint16(J,r),b[V++]=z.getUint16(J+2,r),b[V++]=z.getUint16(J+4,r),J+=m*2}const _=Wt(e,r,u,50714),F=((C=(T=t==null?void 0:t.color_data)==null?void 0:T.dng_levels)==null?void 0:C.dng_cblack)||((k=t==null?void 0:t.color_data)==null?void 0:k.cblack_rawpy_style)||(t==null?void 0:t.black_level_per_channel)||(t==null?void 0:t.cblack),v=Ia(_.length?_:F),w=Wt(e,r,u,50717),P=us(w,(S=(M=t==null?void 0:t.color_data)==null?void 0:M.dng_levels)==null?void 0:S.dng_whitelevel,(R=t==null?void 0:t.color_data)==null?void 0:R.maximum,t==null?void 0:t.white_level)||65535;return{data:b,width:d,height:c,bayerPattern:"RGBG",blackLevels:v,whiteLevel:P,metadata:{...t,format:"DNG_LINEAR_RAW_RGB",description:"Uncompressed DNG LinearRaw RGB",linearRawDngDecoder:!0,bitsPerSample:g.slice(0,m),samplesPerPixel:m,photometric:x},isThreePlane:!0,threePlaneTransfer:"linear"}}async function Na(n){const t=Ra(n);if(!t)return null;let e={};const i=await cs(),r=new i;try{await r.open(n,{}),e=await r.metadata(!0)}catch(s){console.warn("LinearRaw DNG metadata enrichment failed",s)}finally{try{r.delete?r.delete():r.close()}catch{}}try{const s=await ls.parse(n.buffer);s&&(e={...e,...s})}catch(s){console.warn("exifr parsing failed for LinearRaw DNG",s)}return t.metadata={...e,...t.metadata},t}const La=async n=>{var r,s,o,a,l,u,d,c,m,h,f,p,g,x;const t=await Na(n);if(t)return t;const e=await cs(),i=new e;try{if(await i.open(n,{}),typeof i.getRawImage!="function")throw new Error("WASM mismatch");let y={};try{y=await i.metadata(!0)}catch(A){console.warn("Metadata error before raw extraction",A)}const b=i.getRawImage();let _=new Uint16Array(b.data),F={...y};try{const A=await ls.parse(n.buffer);A&&(F={...F,...A})}catch(A){console.warn("exifr parsing failed for RAW buffer",A)}const v=((r=y.idata)==null?void 0:r.filters)||0,w=((s=y.idata)==null?void 0:s.colors)||0,P=v===0&&w===3,T=v===9;let C=[0,0,0,0],k=!1,M=null;if(i.getBlackLevels)try{const A=i.getBlackLevels();M=A,A.dng_cblack&&A.dng_cblack.length===4&&Array.from(A.dng_cblack).some(L=>L>0)?(C=Array.from(A.dng_cblack).map(Number),k=!0):A.cblack&&A.cblack.length===4&&Array.from(A.cblack).some(L=>L>0)?(C=Array.from(A.cblack).map(Number),k=!0):typeof A.black=="number"&&A.black>0&&(C=[A.black,A.black,A.black,A.black],k=!0)}catch(A){console.warn("getBlackLevels binding failed",A)}if(!k){let A=[];if((o=y.color_data)!=null&&o.cblack_rawpy_style)A=y.color_data.cblack_rawpy_style;else if((l=(a=y.color_data)==null?void 0:a.dng_levels)!=null&&l.dng_cblack)A=y.color_data.dng_levels.dng_cblack;else if(((u=y.black_level_per_channel)==null?void 0:u.length)>=4)A=y.black_level_per_channel;else if(((d=y.cblack)==null?void 0:d.length)>=4)A=y.cblack;else if(((m=(c=y.color)==null?void 0:c.cblack)==null?void 0:m.length)>=4)A=y.color.cblack;else{const L=y.black_level||y.color_data&&y.color_data.black||0;A=[L,L,L,L]}C=[Number(A[0])||0,Number(A[1])||0,Number(A[2])||0,Number(A[3])||0]}const S=Number.isFinite(Number(b.bits))&&Number(b.bits)>0?Math.pow(2,Number(b.bits))-1:null,R=us(y.white_level,(f=(h=y.color_data)==null?void 0:h.dng_levels)==null?void 0:f.dng_whitelevel,(p=y.color_data)==null?void 0:p.maximum,M==null?void 0:M.maximum,(g=y.color_data)==null?void 0:g.fmaximum,(x=y.color_data)==null?void 0:x.data_maximum,S)||16383;return{data:_,width:b.width,height:b.height,bayerPattern:y.color_desc||"RGGB",blackLevels:C,whiteLevel:R,metadata:F,isThreePlane:P,threePlaneTransfer:P?"linear":void 0,isXTrans:T}}finally{i.delete?i.delete():i.close()}};async function Ea(n){if(Zs(n)){const i=await sa(n);if(!i)throw new Error("Sony cRAW HQ decoder did not return image data.");return i.rawImageData}const e=new Uint8Array(n);return La(e)}function Ua(n,t,e){const i=hs(n,e),r=Math.floor(t.x),s=Math.floor(t.y),o=Math.floor(t.w),a=Math.floor(t.h),l=new Uint16Array(o*a);for(let u=0;u<a;u++){const d=s+u,c=u*o;for(let m=0;m<o;m++)l[c+m]=i(r+m,d)}return{data:l,width:o,height:a}}function Da(n,t,e){const i=Ba(n,t,e);if(!i)return null;const r=n.width,s=n.height,o=new Uint16Array(r*s);for(let a=0;a<s;a++){const l=a*r;for(let u=0;u<r;u++)o[l+u]=i(u,a)}return{kind:"u16-mono",data:o,width:r,height:s}}function Ba(n,t,e){return t.renderMode==="advanced-zero-dep"&&t.advancedZeroDep?hs(n,t,e):t.renderMode==="zero-dependency"?Oa(n,t,e):null}function hs(n,t,e){if(!t.advancedZeroDep)throw new Error("Unmixing settings not found in DisplaySettings.");const{bg:i,fg:r}=t.advancedZeroDep,s=ds(e,t.advancedZeroDep.bl),{data:o,width:a,whiteLevel:l}=n,u=i.map((f,p)=>Math.max(0,f-s[p])),d=r.map((f,p)=>Math.max(0,f-s[p])),c=(u[1]+u[3])/2,m=(d[1]+d[3])/2,h=Math.pow(2,t.exposure);return(f,p)=>{if(f<0||p<0||f>=a||p>=n.height)return 0;const g=p%2,x=f%2;let y=0;!g&&!x?y=0:!g&&x?y=1:g&&!x?y=3:y=2;const b=o[p*a+f],_=s[y],F=u[y],v=d[y],w=Math.max(b-_,0),P=v-F||1e-9,T=(w-F)/P;let C;return T<0?C=w*(c/Math.max(F,1e-9)):T>1?C=w*(m/Math.max(v,1e-9)):C=(1-T)*c+T*m,C*=h,Math.max(0,Math.min(65535,Math.round(C)))}}function Oa(n,t,e){const{data:i,width:r,height:s,whiteLevel:o}=n,a=ds(e,t.blackLevel||[0,0,0,0]),l=Ga(n.bayerPattern),u=t.wbGains?t.wbGains[0]:1,d=t.wbGains?t.wbGains[1]:1,c=Math.pow(2,t.exposure||0);return(m,h)=>{if(m<0||h<0||m>=r||h>=s)return 0;const f=fs(m,h),p=Xa(l,m,h),g=i[h*r+m],x=a[f];let y=(g-x)/Math.max(1,o-x);return y=Math.max(0,Math.min(1,y)),y*=c,p==="R"?y*=u:p==="B"&&(y*=d),Va(y)}}function ds(n,t){if(typeof n=="number"&&Number.isFinite(n)){const e=Math.max(0,n);return[e,e,e,e]}return Array.isArray(n)&&n.length===4?[Math.max(0,n[0]??0),Math.max(0,n[1]??0),Math.max(0,n[2]??0),Math.max(0,n[3]??0)]:[Math.max(0,t[0]??0),Math.max(0,t[1]??0),Math.max(0,t[2]??0),Math.max(0,t[3]??0)]}function Va(n){return Math.max(0,Math.min(65535,Math.round(Math.max(0,Math.min(1,n))*65535)))}function Ga(n){const t=(n||"RGGB").toUpperCase().trim();return t.length>=4&&/^[RGB]{4}$/.test(t.slice(0,4))?t.slice(0,4):"RGGB"}function fs(n,t){return(t&1)<<1|n&1}function Xa(n,t,e){return n[fs(t,e)]}const ft=(n,t=0)=>({real:n,imag:t}),dn=(n,t)=>({real:n.real+t.real,imag:n.imag+t.imag}),un=(n,t)=>({real:n.real-t.real,imag:n.imag-t.imag}),It=(n,t)=>({real:n.real*t.real-n.imag*t.imag,imag:n.real*t.imag+n.imag*t.real}),Ee=(n,t)=>{const e=t.real*t.real+t.imag*t.imag;return e===0?ft(0):{real:(n.real*t.real+n.imag*t.imag)/e,imag:(n.imag*t.real-n.real*t.imag)/e}},qe=n=>Math.hypot(n.real,n.imag),ps=n=>{const t=qe(n);if(t===0)return ft(0);const e=Math.sqrt(t),i=Math.atan2(n.imag,n.real);return ft(e*Math.cos(i/2),e*Math.sin(i/2))};function za(n,t){const e=n.length-1;if(e<0)return{p:ft(0),dp:ft(0),d2p:ft(0)};let i=ft(n[e].real,n[e].imag),r=ft(0),s=ft(0);for(let o=e-1;o>=0;o--)s=dn(It(r,ft(2)),It(t,s)),r=dn(i,It(t,r)),i=dn(ft(n[o].real,n[o].imag),It(t,i));return{p:i,dp:r,d2p:s}}function Xi(n,t,e=80){const r=n.length-1;if(r<=0)return{root:t,iterations:0};if(r===1)return{root:Ee(It(n[0],ft(-1)),n[1]),iterations:0};let s=ft(t.real,t.imag);for(let o=0;o<e;o++){const{p:a,dp:l,d2p:u}=za(n,s);if(qe(a)<1e-14)return{root:s,iterations:o};const d=Ee(l,a),c=It(d,d),m=un(c,Ee(u,a)),h=ft(r),f=ft(r-1),p=un(It(h,m),It(d,d)),g=ps(It(f,p)),x=dn(d,g),y=un(d,g),b=qe(x)>qe(y)?x:y;if(qe(b)<1e-14)return{root:s,iterations:o};const _=Ee(h,b),F=un(s,_);if(qe(_)<1e-14*qe(F))return{root:F,iterations:o+1};s=F}return{root:s,iterations:e}}function Ya(n,t){const e=n.length-1;if(e<=0)return[ft(0)];if(e===1)return[ft(n[0].real,n[0].imag)];const i=new Array(e);i[e-1]=ft(n[e].real,n[e].imag);for(let r=e-2;r>=0;r--){const s=ft(n[r+1].real,n[r+1].imag),o=i[r+1];i[r]=dn(s,It(t,o))}return i}function Wa(n){const t=n.length-1;if(t<=0)return[];if(t===1)return[Ee(It(n[0],ft(-1)),n[1])];const e=[];let i=n.map(s=>ft(s.real,s.imag)),r=t*5;for(;i.length>2&&r-- >0;){const s=ft(.3+Math.random()*.7,.3+Math.random()*.7),{root:o}=Xi(i,s,100),a=Xi(n,o,20);e.push(a.root);const l=Ya(i,o);if(l.length>=i.length){console.warn("polyDeflate did not reduce degree, breaking");break}i=l}if(i.length===2)e.push(Ee(It(i[0],ft(-1)),i[1]));else if(i.length===3){const s=i[2],o=i[1],a=i[0],l=un(It(o,o),It(It(ft(4),s),a)),u=ps(l),d=It(ft(2),s),c=Ee(un(It(ft(-1),o),u),d),m=Ee(dn(It(ft(-1),o),u),d);e.push(c,m)}return e}function ja(n,t,e){const i=[ft(n),ft(-1),ft(n*t),ft(0),ft(n*e)],r=Wa(i);if(r.length===0)return console.warn("laguerreSmallestPositiveRoot: no roots found"),n;let s=1/0,o=!1;for(const l of r)Math.abs(l.imag)<1e-10&&l.real>0&&l.real<s&&(s=l.real,o=!0);return o?Xi(i,ft(s,0),20).root.real:(console.warn("laguerreSmallestPositiveRoot: no positive real root found"),n)}function Ha(n,t,e){if(Math.abs(t)<1e-10&&Math.abs(e)<1e-10)return n;if(n<1e-10)return 0;if(Math.abs(e)<1e-10){const i=-1/(t*n),r=1/t,s=i*i-4*r;if(s<0)return n;const o=Math.sqrt(s),a=-.5*(i+Math.sign(i)*o),l=a,u=r/a;return l>0&&u>0?Math.min(l,u):l>0?l:u>0?u:n}try{return ja(n,t,e)}catch(i){return console.error("laguerreSmallestPositiveRoot failed:",i),n}}function Qa(n,t,e,i,r){const s=n.x-t.x,o=n.y-t.y,a=Math.hypot(s,o)/Math.max(1e-12,e);if(a<1e-12)return{x:n.x,y:n.y};const l=a*a,u=1+(i+r*l)*l;return{x:s/u+t.x,y:o/u+t.y}}function qa(n,t,e,i,r){const s=n.x-t.x,o=n.y-t.y,a=Math.hypot(s,o)/Math.max(1e-12,e);if(a<1e-12)return{x:t.x,y:t.y};const u=Ha(a,i,r)/a;return{x:t.x+s*u,y:t.y+o*u}}function Ut(n,t){const e=Qa(n,{x:t.principalX,y:t.principalY},t.radiusNorm,t.k1,t.k2);return{x:e.x+(t.correctedOffsetX??0),y:e.y+(t.correctedOffsetY??0)}}function ke(n,t){const e={x:n.x-(t.correctedOffsetX??0),y:n.y-(t.correctedOffsetY??0)};return qa(e,{x:t.principalX,y:t.principalY},t.radiusNorm,t.k1,t.k2)}const Ka=`
attribute vec2 a_position;

void main() {
  gl_Position = vec4(a_position, 0.0, 1.0);
}
`,$a=`
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
`,Ja=`
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
`;class Za{constructor(){bt(this,"canvas",null);bt(this,"gl",null);bt(this,"blurProgram",null);bt(this,"sobelProgram",null);bt(this,"positionBuffer",null);bt(this,"blurUniforms",null);bt(this,"sobelUniforms",null);bt(this,"resources",null);bt(this,"initialized",!1);bt(this,"unavailable",!1);bt(this,"maxTextureSize",0)}compute(t,e,i){if(!this.initialized&&!this.init())return null;const r=this.gl,s=this.blurProgram,o=this.sobelProgram,a=this.blurUniforms,l=this.sobelUniforms;if(!r||!s||!o||!a||!l||!this.positionBuffer||!this.canvas||e<=2||i<=2||e>this.maxTextureSize||i>this.maxTextureSize)return null;const u=this.ensureResources(e,i);if(!u)return null;this.canvas.width=e,this.canvas.height=i,r.viewport(0,0,e,i),r.disable(r.BLEND),r.pixelStorei(r.UNPACK_ALIGNMENT,1),r.pixelStorei(r.PACK_ALIGNMENT,1),r.bindBuffer(r.ARRAY_BUFFER,this.positionBuffer),r.activeTexture(r.TEXTURE0),r.bindTexture(r.TEXTURE_2D,u.sourceTexture),r.texImage2D(r.TEXTURE_2D,0,r.LUMINANCE,e,i,0,r.LUMINANCE,r.UNSIGNED_BYTE,t),r.useProgram(s),r.enableVertexAttribArray(0),r.vertexAttribPointer(0,2,r.FLOAT,!1,0,0),r.uniform2f(a.size,e,i),r.uniform1i(a.source,0),r.bindFramebuffer(r.FRAMEBUFFER,u.blurFramebuffer),r.bindTexture(r.TEXTURE_2D,u.sourceTexture),r.drawArrays(r.TRIANGLES,0,6);const d=new Uint8Array(e*i*4);r.readPixels(0,0,e,i,r.RGBA,r.UNSIGNED_BYTE,d),r.useProgram(o),r.uniform2f(l.size,e,i),r.uniform1i(l.blurred,0),r.bindFramebuffer(r.FRAMEBUFFER,u.sobelFramebuffer),r.bindTexture(r.TEXTURE_2D,u.blurTexture),r.drawArrays(r.TRIANGLES,0,6);const c=new Uint8Array(e*i*4);r.readPixels(0,0,e,i,r.RGBA,r.UNSIGNED_BYTE,c),r.disableVertexAttribArray(0),r.bindFramebuffer(r.FRAMEBUFFER,null),r.bindBuffer(r.ARRAY_BUFFER,null),r.bindTexture(r.TEXTURE_2D,null);const m=new Uint8Array(e*i);for(let g=0,x=0;g<m.length;g++,x+=4)m[g]=d[x];const h=new Float32Array(e*i),f=new Float32Array(e*i),p=new Float32Array(e*i);for(let g=0,x=0;g<h.length;g++,x+=4){const y=(c[x]|c[x+1]<<8)-32768,b=(c[x+2]|c[x+3]<<8)-32768;h[g]=y,f[g]=b,p[g]=Math.sqrt(y*y+b*b)}return{blurredGray:m,gx:h,gy:f,magnitude:p}}init(){if(this.initialized&&this.gl&&this.blurProgram&&this.sobelProgram)return!0;if(this.unavailable)return!1;const t=this.createCanvas();if(!t)return this.unavailable=!0,!1;const e=t.getContext("webgl",{alpha:!1,antialias:!1,depth:!1,stencil:!1,premultipliedAlpha:!1,preserveDrawingBuffer:!1});if(!e)return this.unavailable=!0,!1;const i=this.compileShader(e,e.VERTEX_SHADER,Ka),r=this.compileShader(e,e.FRAGMENT_SHADER,$a),s=this.compileShader(e,e.FRAGMENT_SHADER,Ja);if(!i||!r||!s)return i&&e.deleteShader(i),r&&e.deleteShader(r),s&&e.deleteShader(s),this.unavailable=!0,!1;const o=this.createProgram(e,i,r),a=this.createProgram(e,i,s);if(e.deleteShader(i),e.deleteShader(r),e.deleteShader(s),!o||!a)return o&&e.deleteProgram(o),a&&e.deleteProgram(a),this.unavailable=!0,!1;const l=e.createBuffer();return l?(e.bindBuffer(e.ARRAY_BUFFER,l),e.bufferData(e.ARRAY_BUFFER,new Float32Array([-1,-1,1,-1,-1,1,-1,1,1,-1,1,1]),e.STATIC_DRAW),e.bindBuffer(e.ARRAY_BUFFER,null),this.canvas=t,this.gl=e,this.blurProgram=o,this.sobelProgram=a,this.positionBuffer=l,this.blurUniforms={source:e.getUniformLocation(o,"u_source"),size:e.getUniformLocation(o,"u_size")},this.sobelUniforms={blurred:e.getUniformLocation(a,"u_blurred"),size:e.getUniformLocation(a,"u_size")},this.maxTextureSize=Number(e.getParameter(e.MAX_TEXTURE_SIZE)||0),this.initialized=!0,!0):(e.deleteProgram(o),e.deleteProgram(a),this.unavailable=!0,!1)}createCanvas(){return typeof OffscreenCanvas<"u"?new OffscreenCanvas(1,1):typeof document<"u"?document.createElement("canvas"):null}ensureResources(t,e){const i=this.gl;if(!i)return null;if(this.resources&&this.resources.width===t&&this.resources.height===e)return this.resources;this.disposeResources();const r=this.createTexture(i.LUMINANCE,t,e,i.LUMINANCE,i.UNSIGNED_BYTE,null),s=this.createTexture(i.RGBA,t,e,i.RGBA,i.UNSIGNED_BYTE,null),o=this.createTexture(i.RGBA,t,e,i.RGBA,i.UNSIGNED_BYTE,null),a=this.createFramebuffer(s),l=this.createFramebuffer(o);return!r||!s||!o||!a||!l?(r&&i.deleteTexture(r),s&&i.deleteTexture(s),o&&i.deleteTexture(o),a&&i.deleteFramebuffer(a),l&&i.deleteFramebuffer(l),null):(this.resources={width:t,height:e,sourceTexture:r,blurTexture:s,sobelTexture:o,blurFramebuffer:a,sobelFramebuffer:l},this.resources)}createTexture(t,e,i,r,s,o){const a=this.gl;if(!a)return null;const l=a.createTexture();return l?(a.bindTexture(a.TEXTURE_2D,l),a.texParameteri(a.TEXTURE_2D,a.TEXTURE_MIN_FILTER,a.NEAREST),a.texParameteri(a.TEXTURE_2D,a.TEXTURE_MAG_FILTER,a.NEAREST),a.texParameteri(a.TEXTURE_2D,a.TEXTURE_WRAP_S,a.CLAMP_TO_EDGE),a.texParameteri(a.TEXTURE_2D,a.TEXTURE_WRAP_T,a.CLAMP_TO_EDGE),a.texImage2D(a.TEXTURE_2D,0,t,e,i,0,r,s,o),a.bindTexture(a.TEXTURE_2D,null),l):null}createFramebuffer(t){const e=this.gl;if(!e||!t)return null;const i=e.createFramebuffer();if(!i)return null;e.bindFramebuffer(e.FRAMEBUFFER,i),e.framebufferTexture2D(e.FRAMEBUFFER,e.COLOR_ATTACHMENT0,e.TEXTURE_2D,t,0);const r=e.checkFramebufferStatus(e.FRAMEBUFFER);return e.bindFramebuffer(e.FRAMEBUFFER,null),r!==e.FRAMEBUFFER_COMPLETE?(e.deleteFramebuffer(i),null):i}compileShader(t,e,i){const r=t.createShader(e);return r?(t.shaderSource(r,i),t.compileShader(r),t.getShaderParameter(r,t.COMPILE_STATUS)?r:(console.error("[SFR Auto Detect WebGL] shader compile failed",t.getShaderInfoLog(r)),t.deleteShader(r),null)):null}createProgram(t,e,i){const r=t.createProgram();return r?(t.attachShader(r,e),t.attachShader(r,i),t.bindAttribLocation(r,0,"a_position"),t.linkProgram(r),t.getProgramParameter(r,t.LINK_STATUS)?r:(console.error("[SFR Auto Detect WebGL] program link failed",t.getProgramInfoLog(r)),t.deleteProgram(r),null)):null}disposeResources(){const t=this.gl,e=this.resources;if(!t||!e){this.resources=null;return}t.deleteTexture(e.sourceTexture),t.deleteTexture(e.blurTexture),t.deleteTexture(e.sobelTexture),t.deleteFramebuffer(e.blurFramebuffer),t.deleteFramebuffer(e.sobelFramebuffer),this.resources=null}}const to=new Za,zi={gradientPercentiles:[.82,.88,.92,.95,.98,.995],downsampleMaxSide:1600,minComponentAreaRatio:15e-6,maxComponentAreaRatio:.35,minComponentAreaPx:20,minEdgePoints:24,extentQuantileLow:.02,extentQuantileHigh:.98,cornerTrimRatio:.18,minSpanPx:8,maxAspectRatio:2,bandScale:.16,bandMinPx:1.75,bandMaxPx:14,minPointContrast:6,minSidePoints:3,minCoverageRatio:.15,minCenterCoverageRatio:.2,filterBlockPurity:!0,innerPurityStdScale:1.5,outerMeanSpreadLimit:51,minAxisDot:.6,residualLimitFloor:.01,residualLimitScale:.25,minQuadArea:48,minSideLength:10,minOuterContrast:5,sampleHalfWidthRatio:.25};function eo(n,t,e,i,r,s){const o=n.width,a=n.height,l=(n.bayerPattern||"RGGB").toUpperCase(),u=[],d=[],c=s!=null&&s.correctedRect?Lt*2:Lt,m=Math.max(1,Math.min(r,c)),h=e.p2.x-e.p1.x,f=e.p2.y-e.p1.y,p=Math.hypot(h,f);if(!Number.isFinite(p)||p<=1e-6)return null;const g=h/p,x=f/p,y=-x,b=g,_={x:(e.p1.x+e.p2.x)*.5,y:(e.p1.y+e.p2.y)*.5},F=s!=null&&s.correctedRect?Mt(s.correctedRect,o,a):Mt(Vt(be(e,c*4+2)??[e.p1,e.p2],2),o,a);if(!F)return null;const v=(s==null?void 0:s.correctedScanlinesOverride)??(s!=null&&s.distortedRect?rr(Mt(s.distortedRect,o,a)??s.distortedRect,t,o,a):Ts(F,e,Math.max(1,i),m*4+.5,o,a));if(!v||v.size===0)return null;const w=Is(v,t,o,a);if(w.size===0)return null;const P=!sr(t);for(const[C,k]of w){if(C<0||C>=a)continue;const M=C*o;for(let S=k.start;S<=k.end;S++){if(S<0||S>=o||!vt(S,C,l,s==null?void 0:s.greenPhase))continue;const R={x:S,y:C},A=Ut(R,t);if(!Number.isFinite(A.x)||!Number.isFinite(A.y)||Math.round(A.x)<0||Math.round(A.x)>=o||Math.round(A.y)<0||Math.round(A.y)>=a)continue;const L=A.x-_.x,X=A.y-_.y,D=L*g+X*x;let E=L*y+X*b;if(P){const U=Rs(D,g,x,_,R,t);if(!U)continue;const O=.5*(U.a+U.b),z=Kt(U.a,g,x,_,t),J=Kt(O,g,x,_,t),V=Kt(U.b,g,x,_,t),Q=or({x:U.a,y:Math.hypot(z.x-R.x,z.y-R.y)},{x:O,y:Math.hypot(J.x-R.x,J.y-R.y)},{x:U.b,y:Math.hypot(V.x-R.x,V.y-R.y)});if(!Number.isFinite(Q))continue;const W=ar(Q,g,x,_,t),tt=Math.hypot(W.x,W.y);if(!Number.isFinite(tt)||tt<=1e-9)continue;const at=W.x/tt,Z=-(W.y/tt),q=at,I=Kt(Q,g,x,_,t);E=(R.x-I.x)*Z+(R.y-I.y)*q}!Number.isFinite(D)||Math.abs(D)>Math.max(1,i)||!Number.isFinite(E)||Math.abs(E)>m||(u.push(E),d.push(Math.max(0,n.data[M+S]-yn(s==null?void 0:s.blackLevel,S,C))))}}if(u.length<8)return null;const T=Math.abs(h)>=Math.abs(f)?1:2;return s!=null&&s.forceLegacyModel?Fn(u,d,T,c):xn(u,d,T,c)}function Zi(n,t){const e=n.length;let i=0,r=0,s=0,o=0;for(let l=0;l<e;l++)i+=n[l],r+=t[l],s+=n[l]*t[l],o+=n[l]*n[l];const a=e*o-i*i;return a===0?{slope:0,intercept:0}:{slope:(e*s-i*r)/a,intercept:(r*o-i*s)/a}}function no(n,t){const e=n.length,i=new Array(e).fill(0),r=2*t;for(let s=0;s<e;s++){const o=s>0?n[s-1]:n[0],a=s<e-1?n[s+1]:n[e-1];i[s]=(a-o)/r}return i}const te=-1e7,fn=13,ae=512,Et=8,ci=1/Et,Lt=28,io=[[0,0,0,0,0,-.085714285714286,.342857142857143,.485714285714286,.342857142857143,-.085714285714286,0,0,0,0,0],[0,0,0,0,-.095238095238095,.142857142857143,.285714285714286,.333333333333333,.285714285714286,.142857142857143,-.095238095238095,0,0,0,0],[0,0,0,-.090909090909091,.060606060606061,.168831168831169,.233766233766234,.255411255411255,.233766233766234,.168831168831169,.060606060606061,-.090909090909091,0,0,0],[0,0,-.083916083916084,.020979020979021,.102564102564103,.160839160839161,.195804195804196,.207459207459208,.195804195804196,.160839160839161,.102564102564103,.020979020979021,-.083916083916084,0,0],[0,-.076923076923077,0,.062937062937063,.111888111888112,.146853146853147,.167832167832168,.174825174825175,.167832167832168,.146853146853147,.111888111888112,.062937062937063,0,-.076923076923077,0],[-.070588235294118,-.011764705882353,.038009049773756,.078733031674208,.110407239819004,.133031674208145,.146606334841629,.151131221719457,.146606334841629,.133031674208145,.110407239819004,.078733031674208,.038009049773756,-.011764705882353,-.070588235294118]];function pn(n,t,e,i=1){const r=Math.max(1e-6,e*.5),s=Math.max(1e-6,t*i),o=Math.exp(-s*r),a=1-o;if(!Number.isFinite(a)||Math.abs(a)<=1e-9)return Math.abs(n)<=r?1:0;if(Math.abs(n)<r){const l=2-2*o*Math.cosh(s*n),u=2*Math.sinh(s*r)*a;return!Number.isFinite(u)||Math.abs(u)<=1e-9?0:l/u}return Math.exp(-s*Math.abs(n))/a}function ms(n,t,e,i,r,s){const o=n.length;if(o===0)return[];if(s<1)return n;const a=Math.min(s,32),l=new Array(o).fill(0);l[0]=n[0];for(let h=1;h<o;h++)l[h]=l[h-1]+n[h];const u=(h,f)=>{const p=Math.max(0,h),g=Math.min(o-1,f);return g<p?n[Math.max(0,Math.min(o-1,h))]??0:(l[g]-(p>0?l[p-1]:0))/(g-p+1)},d=a*2,c=a,m=1;for(let h=Math.max(t+c,i-d);h<i;h++){const f=Math.max(m,Math.trunc((i-h)*c/Math.max(1,d)));n[h]=u(h-f,h+f)}for(let h=Math.min(r+d-1,e-c-1);h>r;h--){const f=Math.max(m,Math.trunc((h-r)*c/Math.max(1,d)));n[h]=u(h-f,h+f)}for(let h=c+1;h<i-d;h++)n[h]=u(h-c,h+c);for(let h=Math.min(r+d,e-c-1);h<o-c-1;h++)n[h]=u(h-c,h+c);return n}function gs(n){return!Number.isFinite(n)||Math.abs(n)<=1e-9?1:Math.sin(n)/n}let zn=null,Ei=null;function ro(){if(zn)return zn;const n=.625,t=1/128,e=Math.max(16,Math.round(n*2/t)+1),i=[],r=[];for(let c=0;c<e;c++){const m=-n+c*t;i.push(m),r.push(Math.abs(m)<=n?pn(m,fn,.125,1):0)}const s=4,o=1/1024,a=Math.round(s/o)+1,l=new Array(a).fill(0),u=new Array(a).fill(1);let d=0;for(let c=0;c<r.length;c++)d+=r[c];d=Math.max(1e-9,d);for(let c=0;c<a;c++){const m=c*o;l[c]=m;let h=0;for(let f=0;f<i.length;f++)h+=r[f]*Math.cos(2*Math.PI*m*i[f]);u[c]=Math.max(1e-6,Math.abs(h)/d)}return zn={freqs:l,values:u},zn}function so(n,t){const e=Math.max(1e-6,gs(Math.PI*n*t)),i=ro(),r=Ze(Math.max(0,Math.min(i.freqs[i.freqs.length-1],n)),i.freqs,i.values);return Math.max(1e-6,e*r)}function ao(){if(Ei)return Ei;const n=new Array(ae/16*4).fill(1),t=ae*16,e=new Float32Array(t);for(let s=0;s<t;s++){const o=(s-t/2)/(16*Et);e[s]=Math.abs(o)<=.625?pn(o,fn,ci,1):0}const i=new kn(t);i.transform(e);const r=Math.max(1e-9,Math.abs(i._real[0]));n[0]=1;for(let s=1;s<n.length;s++){const o=gs(Math.PI*s/256),a=Math.max(1e-6,Math.hypot(i._real[s],i._imag[s])/r);n[s]=Math.max(1e-6,o*a)}return Ei=n,n}function oo(n,t,e){if(n.length===0||t.length!==n.length||!(e>0))return null;const i=n[0],r=n[n.length-1],s=Math.floor((r-i)/e);if(s<16)return null;const o=Math.max(0,Math.min(s-1,Math.round(-i/e))),a=Math.max(2,Math.round(2/e)),l=Math.max(1,Math.round(.5/e)),u=5,d=.125/Math.max(e,1e-6),c=B=>i+B*e,h=((B,Z)=>{let q=Math.max(0,B),I=Math.min(s-1,Z);if(I-q<8)return null;let H=0;for(;;){const G=new Array(s).fill(0),Y=new Array(s).fill(0);for(let $=0;$<n.length;$++){const Ft=Math.max(0,Math.min(s-1,Math.trunc((n[$]-i)/e))),Tt=Math.max(q,Ft-u),xt=Math.min(I-1,Ft+u);for(let wt=Tt;wt<=xt;wt++){const ne=c(wt),Jt=Math.max(0,1-Math.abs((n[$]-ne)*1.75*d));Jt<=0||(G[wt]+=t[$]*Jt,Y[wt]+=Jt)}}const et=new Array(s).fill(te);let nt=0,it=0,ct=0,dt=0,rt=-1,Dt=-1;const Gt=Math.max(o-Math.round(s/8),q+a),j=Math.min(o+Math.round(s/8),I-a);for(let $=Math.max(0,q-1);$<=Math.min(s-1,I+1);$++)Y[$]>0&&(et[$]=G[$]/Y[$],$<Gt&&(nt+=et[$],ct++),$>j&&(it+=et[$],dt++),rt<0&&(rt=$),Dt=$);if(rt<0||Dt<0||ct===0||dt===0)return null;for(let $=rt-1;$>=0;$--)et[$]=et[rt];for(let $=Dt+1;$<s;$++)et[$]=et[Dt];const ut=Math.max(2,a),ot=new Array(s).fill(0);let pt=o;for(let $=ut+1;$<s-1-ut;$++){let Ft=0,Tt=0;for(let xt=-ut;xt<=ut;xt++)Tt+=et[$+xt]*xt,Ft+=xt*xt;ot[$]=Ft>0?Tt/Ft:0,Math.abs(ot[$])>Math.abs(ot[pt]??0)&&$>q+ut&&$<I-ut-1&&(pt=$)}const mt=Math.max(1,Math.round(2/e)),ee=Math.max(mt+1,Math.round(12/e)),Bt=Math.abs(pt-o);if(Bt>mt&&Bt<ee)return null;let _t=0;for(let $=Math.max(0,o-ut);$<=Math.min(s-1,o+ut);$++)Math.abs(ot[$])>Math.abs(_t)&&(_t=ot[$]);if(!Number.isFinite(_t)||Math.abs(_t)<=1e-9)return null;const ce=Math.abs(_t*.001);let $t=!1,Fe=!1;for(let $=Math.max(0,o-a);$>=q+l;$--)if(ot[$]*_t<0&&Math.abs(ot[$])>ce){let Ft=0,Tt=0,xt=0;for(let wt=$;wt>=q;wt--)ot[wt]*_t<0&&(Ft++,Tt=Math.max(Tt,Math.abs(ot[wt]))),xt++;if(Ft>xt*.4&&Tt/Math.abs(_t)>.25||Ft>.9*xt&&xt>a){q=Math.min($,o-a),$t=!0;break}}for(let $=Math.min(s-1,o+a);$<I-l;$++)if(ot[$]*_t<0&&Math.abs(ot[$])>ce){let Ft=0,Tt=0,xt=0;for(let wt=$;wt<I;wt++)ot[wt]*_t<0&&(Ft++,Tt=Math.max(Tt,Math.abs(ot[wt]))),xt++;if(Ft>xt*.4&&Tt/Math.abs(_t)>.25||Ft>.9*xt&&xt>a){I=Math.max($,o+a),$t=!0;break}}if($t&&((o-q<Math.max(1,Math.round(4/e))||I-o<Math.max(1,Math.round(4/e)))&&(Fe=!0),H<2)){H++;continue}return{sampled:et,fftLeft:q,fftRight:I,leftMean:nt/ct,rightMean:it/dt,peakSlopeIdx:pt,slopes:ot,clipped:$t,dodgy:Fe}}})(0,s-1);if(!h)return null;const f=2,p=Math.max(h.leftMean,h.rightMean),g=Math.min(h.leftMean,h.rightMean);let x=h.fftLeft,y=h.fftLeft,b=1/0,_=1/0;for(let B=h.fftLeft;B<=h.fftRight;B++){let Z=0,q=0;for(let Y=-f;Y<=f;Y++){const et=h.sampled[Math.max(0,Math.min(s-1,B+Y))];Z+=et,q++}const I=Z/Math.max(1,q),H=Math.abs(I-g-.1*(p-g)),G=Math.abs(I-g-.9*(p-g));H<b&&(b=H,x=B),G<_&&(_=G,y=B)}const F=Math.max(4,Math.abs(x-y)*e);if(x<y){const B=x;x=y,y=B}const v=Math.max(l,l+2*Math.round(F/Math.max(e,1e-6)));x+=v,y-=v;const w=Math.max(Math.abs(x-o),Math.abs(y-o),Math.max(a,Math.round(4/Math.max(e,1e-6)))),P=1.85,T=.5,C=new Array(s).fill(0),k=new Array(s).fill(0);for(let B=0;B<n.length;B++){const Z=Math.max(0,Math.min(s-1,Math.trunc((n[B]-i)/e)));let q=5;Math.abs(Z-o)>T*w&&(q=Math.abs(Z-o)>2*T*w?12:7);const I=Math.max(h.fftLeft,Z-q),H=Math.min(h.fftRight-1,Z+q);if(H<o-P*w||I>o+P*w){for(let G=I;G<=H;G++)C[G]+=t[B],k[G]+=1;continue}for(let G=I;G<=H;G++){let Y=1;if(Math.abs(G-o)<P*w){const et=c(G);if(Math.abs(G-o)<w*T)Y=pn(n[B]-et,fn,e,1);else{const nt=(Math.abs(G-o)/Math.max(1e-6,w)-T)/Math.max(1e-6,P-T),it=1*(1-nt)+.01*nt;Y=pn(n[B]-et,fn,e,it)}}!(Y>0)||!Number.isFinite(Y)||(C[G]+=t[B]*Y,k[G]+=Y)}}const M=new Array(s).fill(0);let S=-1,R=-1;for(let B=0;B<s;B++)k[B]>0?(M[B]=C[B]/k[B],S<0&&(S=B),R=B):M[B]=te;if(S<0||R<0)return null;const A=Math.max(1,Math.round(3/Math.max(e,1e-6)));let L=M[S],X=1;for(let B=S+1;B<o&&X<A;B++)M[B]!==te&&(L+=M[B],X++);L/=X;let D=M[R],E=1;for(let B=R-1;B>o&&E<A;B--)M[B]!==te&&(D+=M[B],E++);D/=E;for(let B=S-1;B>=0;B--)M[B]=L;for(let B=R+1;B<s;B++)M[B]=D;for(let B=S+1;B<R;B++){if(M[B]!==te)continue;let Z=B-1;for(;Z>=0&&M[Z]===te;)Z--;let q=B+1;for(;q<s&&M[q]===te;)q++;if(Z>=0&&q<s){const I=(B-Z)/Math.max(1,q-Z);M[B]=M[Z]*(1-I)+M[q]*I}else Z>=0?M[B]=M[Z]:q<s&&(M[B]=M[q])}const O=[...L<=D?M:[...M].reverse()],z=L<=D?O:O.reverse(),J=Math.max(Math.round(o-P*w),h.fftLeft+2),V=Math.min(Math.round(o+P*w),h.fftRight-3),Q=Math.max(1,Math.round(2/Math.max(e,1e-6))),W=ms(z,h.fftLeft,h.fftRight,J,V,Q),tt=new Array(s).fill(0);let at=W[Math.max(0,Math.min(s-1,h.fftLeft))]??W[0]??0;for(let B=h.fftLeft;B<=h.fftRight;B++){const Z=W[B]??at,q=W[Math.min(s-1,B+1)]??Z;tt[B]=q-at,at=Z}return{esf:W,lsfFull:tt}}function ys(n,t,e=Lt){if(n.length===0||t.length!==n.length)return null;const i=ae,r=i/2,s=ci,o=2*Et,a=Math.max(1,Math.round(.5*Et)),l=5,u=Math.max(0,Math.round(r-e*Et)),d=Math.min(i-1,Math.round(r+e*Et));if(d-u<32)return null;let c=new Array(i).fill(te),m=0,h=0,f=0,p=0,g=u,x=d,y=u,b=d,_=0;for(;;){const I=new Array(i).fill(0),H=new Array(i).fill(0);c=new Array(i).fill(te),m=0,h=0,f=0,p=0;let G=-1,Y=-1;for(let j=0;j<n.length;j++){const ut=Math.trunc(n[j]*Et+r),ot=Math.max(y,ut-l),pt=Math.min(b-1,ut+l);for(let mt=ot;mt<=pt;mt++){const ee=(mt-r)*s,Bt=Math.max(0,1-Math.abs((n[j]-ee)*1.75));Bt>0&&(H[mt]+=t[j]*Bt,I[mt]+=Bt)}}const et=Math.max(r-i/8,y+2*Et),nt=Math.min(r+i/8,b-2*Et);for(let j=Math.max(0,y-1);j<=Math.min(i-1,b+1);j++)I[j]>0&&(c[j]=H[j]/I[j],j<et&&(m+=c[j],f++),j>nt&&(h+=c[j],p++),G<0&&(G=j),Y=j);if(G<0||Y<0||f<=0||p<=0)return null;for(let j=G-1;j>=0;j--)c[j]=c[G];for(let j=Y+1;j<i;j++)c[j]=c[Y];const it=new Array(i).fill(0);let ct=r;const dt=2*Et;for(let j=dt+1;j<i-1-dt;j++){let ut=0,ot=0;for(let pt=-dt;pt<=dt;pt++)ot+=c[j+pt]*pt,ut+=pt*pt;it[j]=ut>0?ot/ut:0,Math.abs(it[j])>Math.abs(it[ct]??0)&&j>y+dt&&j<b-dt-1&&(ct=j)}if(Math.abs(ct-r)>2*Et&&Math.abs(ct-r)<12*Et)return null;let rt=0;for(let j=Math.max(0,r-dt);j<=Math.min(i-1,r+dt);j++)Math.abs(it[j])>Math.abs(rt)&&(rt=it[j]);if(!Number.isFinite(rt)||Math.abs(rt)<=1e-9)return null;const Dt=Math.abs(rt*.001);g=y,x=b;let Gt=!1;for(let j=r-o;j>=y+a;j--)if(it[j]*rt<0&&Math.abs(it[j])>Dt){let ut=0,ot=0,pt=0;for(let mt=j;mt>=y;mt--)it[mt]*rt<0&&(ut++,ot=Math.max(ot,Math.abs(it[mt]))),pt++;if(ut>pt*.4&&ot/Math.abs(rt)>.25||ut>.9*pt&&pt>o){g=Math.min(j,r-o),Gt=!0;break}}for(let j=r+o;j<b-a;j++)if(it[j]*rt<0&&Math.abs(it[j])>Dt){let ut=0,ot=0,pt=0;for(let mt=j;mt<b;mt++)it[mt]*rt<0&&(ut++,ot=Math.max(ot,Math.abs(it[mt]))),pt++;if(ut>pt*.4&&ot/Math.abs(rt)>.25||ut>.9*pt&&pt>o){x=Math.max(j,r+o),Gt=!0;break}}if(Gt&&_<2){y=g,b=x,_++;continue}break}const F=Math.max(m/f,h/p),v=Math.min(m/f,h/p);let w=g,P=g,T=1/0,C=1/0;for(let I=g;I<=x;I++){const H=(c[Math.max(0,I-2)]+c[Math.max(0,I-1)]+c[I]+c[Math.min(i-1,I+1)]+c[Math.min(i-1,I+2)])/5,G=Math.abs(H-v-.1*(F-v)),Y=Math.abs(H-v-.9*(F-v));G<T&&(T=G,w=I),Y<C&&(C=Y,P=I)}if(w<P){const I=w;w=P,P=I}const k=Math.max(4,Math.abs(w-P)*s),M=Math.max(a,a+2*Math.trunc(k/Math.max(s,1e-6)));w+=M,P-=M;const S=Math.max(Math.abs(w-r),Math.abs(P-r),Math.max(o,Math.trunc(4/Math.max(s,1e-6)))),R=new Array(i).fill(0),A=new Array(i).fill(0),L=1.85,X=.5;for(let I=0;I<n.length;I++){const H=Math.trunc(n[I]*Et+r);let G=5;Math.abs(H-r)>X*S&&(G=Math.abs(H-r)>2*X*S?12:7);const Y=Math.max(g,H-G),et=Math.min(x-1,H+G);if(et<r-L*S||Y>r+L*S){for(let nt=Y;nt<=et;nt++)R[nt]+=t[I],A[nt]+=1;continue}for(let nt=Y;nt<=et;nt++){let it=1;if(Math.abs(nt-r)<L*S){const ct=(nt-r)*s;if(Math.abs(nt-r)<S*X)it=pn(n[I]-ct,fn,s,1);else{const dt=(Math.abs(nt-r)/Math.max(1e-6,S)-X)/Math.max(1e-6,L-X),rt=1*(1-dt)+.01*dt;it=pn(n[I]-ct,fn,s,rt)}}!(it>0)||!Number.isFinite(it)||(R[nt]+=t[I]*it,A[nt]+=it)}}const D=new Array(i).fill(0);let E=-1,U=-1;for(let I=Math.max(0,g-1);I<=Math.min(i-1,x+1);I++)A[I]>0?(D[I]=R[I]/A[I],E<0&&(E=I),U=I):D[I]=te;if(E<0||U<0)return null;const O=3*Et;let z=D[E],J=1;for(let I=E+1;I<r&&J<O;I++)D[I]!==te&&(z+=D[I],J++);z/=Math.max(1,J);let V=D[U],Q=1;for(let I=U-1;I>r&&Q<O;I--)D[I]!==te&&(V+=D[I],Q++);V/=Math.max(1,Q);for(let I=E-1;I>=0;I--)D[I]=z;for(let I=U+1;I<i;I++)D[I]=V;const W=Math.max(Math.trunc(r-L*S),g+2),tt=Math.min(Math.trunc(r+L*S),x-3),at=Math.max(1,Math.trunc(2/Math.max(s,1e-6))),B=ms(D,g,x,W,tt,at),Z=new Array(i).fill(0);let q=B[Math.max(0,Math.min(i-1,g))]??B[0]??0;for(let I=g;I<=x;I++){const H=B[I]??q;Z[I]=(B[Math.min(i-1,I+1)]??H)-q,q=H}return{esf:B,lsfFull:Z}}function lo(n){const t=new Array(n.length).fill(0);if(n.length===0)return t;t[0]=n[0];for(let e=1;e<n.length;e++){const i=tr(n[e]-n[e-1]);t[e]=t[e-1]+i}return t}function tr(n){if(!Number.isFinite(n))return n;let t=(n+Math.PI)%(2*Math.PI);return t<0&&(t+=2*Math.PI),t-Math.PI}function co(n,t,e=0){if(n.length===0)return[];const i=Number.isFinite(e)?e:0,r=n.map((o,a)=>{const l=t[a]??0,u=-2*Math.PI*i*l;return tr(o-u)});return lo(r).map((o,a)=>{const l=t[a]??0,u=-2*Math.PI*i*l;return o+u})}function uo(n,t,e,i=.05,r=Number.POSITIVE_INFINITY){const s=Math.min(n.length,t.length);if(s<2)return null;let o=0;if(e)for(let p=1;p<s;p++){const g=e[p];Number.isFinite(g)&&g>o&&(o=g)}const a=e&&o>0?Math.max(1e-6,o*i):0;let l=0,u=0,d=0,c=0,m=0,h=0;for(let p=1;p<s;p++){const g=n[p],x=t[p];if(!Number.isFinite(g)||!Number.isFinite(x)||Math.abs(g)<=1e-12||g>r)continue;const y=e?e[p]:1;if(!Number.isFinite(y)||y<=a)continue;const b=e?y*y:1;l+=b,u+=b*g,d+=b*x,c+=b*g*g,m+=b*g*x,h++}if(h<4||l<=0)return null;const f=l*c-u*u;return Math.abs(f)<=1e-12?null:{slope:(l*m-u*d)/f,intercept:(d*c-u*m)/f,used:h,threshold:a}}function ho(n,t,e=Number.POSITIVE_INFINITY){const i=[],r=[],s=Math.min(n.length,t.length);for(let o=1;o<s;o++)Number.isFinite(n[o])&&Number.isFinite(t[o])&&Math.abs(n[o])>1e-12&&n[o]<=e&&(i.push(n[o]),r.push(t[o]));return i.length<2?{slope:0,intercept:Number.isFinite(t[0])?t[0]:0,used:i.length}:{...Zi(i,r),used:i.length}}function fo(n,t,e,i=.05,r=Number.POSITIVE_INFINITY,s=0){const o=Math.min(n.length,t.length);if(o<4)return null;let a=0;if(e)for(let b=1;b<o;b++){const _=e[b];Number.isFinite(_)&&_>a&&(a=_)}const l=e&&a>0?Math.max(1e-6,a*i):0,u=[];for(let b=1;b<o;b++){const _=n[b],F=t[b];if(!Number.isFinite(_)||!Number.isFinite(F)||Math.abs(_)<=1e-12||_>r)continue;const v=e?e[b]:1;!Number.isFinite(v)||v<=l||u.push({freq:_,phase:F,weight:e?v*v:1})}if(u.length<4)return null;const d=b=>{let _=0,F=0;for(const C of u){const k=C.phase+2*Math.PI*b*C.freq;_+=C.weight*Math.sin(k),F+=C.weight*Math.cos(k)}const v=Math.atan2(_,F);let w=0,P=0;const T=.65;for(const C of u){const k=C.phase+2*Math.PI*b*C.freq,M=Math.abs(tr(k-v)),S=M<=T?M*M:T*(2*M-T);w+=C.weight*S,P+=C.weight}return{score:P>0?w/P:Number.POSITIVE_INFINITY,intercept:v}},c=Number.isFinite(s)?s:0,m=Math.max(2,Math.min(8,Math.abs(c)>1e-6?4:2)),h=.02;let f=c,p=d(f);for(let b=c-m;b<=c+m+h*.5;b+=h){const _=d(b);_.score<p.score&&(p=_,f=b)}let g=f-h*2,x=f+h*2;for(let b=0;b<32;b++){const _=g+(x-g)/3,F=x-(x-g)/3,v=d(_).score,w=d(F).score;v<w?x=F:g=_}const y=(g+x)*.5;return p=d(y),{slope:-2*Math.PI*y,intercept:p.intercept,used:u.length,threshold:l}}function po(n,t){if(Number.isFinite(t)&&t>0)return t;let e=0;for(const i of n)Number.isFinite(i)&&i>e&&(e=i);return e>0?e:Number.POSITIVE_INFINITY}function xs(n,t,e,i=Number.POSITIVE_INFINITY,r=0){const s=co(n,t,r),o=po(t,i),a=fo(t,n,e,.05,o,r),l=a?null:uo(t,s,e,.05,o),u=a||l?null:ho(t,s,o),d=(a==null?void 0:a.slope)??(l==null?void 0:l.slope)??(u==null?void 0:u.slope)??0,c=(a==null?void 0:a.intercept)??(l==null?void 0:l.intercept)??(u==null?void 0:u.intercept)??0,m=s.map((y,b)=>y-(d*(t[b]??0)+c)),h=Number.isFinite(m[0])?m[0]:0,f=t.map(y=>d*y+c+h),p=m.map(y=>y-h),g=c+h,x=Number.isFinite(d)?-d/(2*Math.PI):null;return{raw:[...n],unwrapped:s,linear:f,residual:p,fit:{groupDelayPx:x===null?null:x-r,absoluteGroupDelayPx:x,referenceDelayPx:r,slopeRadPerCycle:Number.isFinite(d)?d:null,interceptRad:Number.isFinite(g)?g:null,fitPointCount:(a==null?void 0:a.used)??(l==null?void 0:l.used)??(u==null?void 0:u.used)??0,fitWeightThreshold:(a==null?void 0:a.threshold)??(l==null?void 0:l.threshold)??0,fitDomain:"cycles-per-pixel",fitMaxFreqCyclesPerPixel:o}}}function bs(n,t,e){const i=[],r=[],s=[],o=[];for(let a=0;a<e.length;a++){const l=e[a];i.push(Ze(l,t,n.raw)),r.push(Ze(l,t,n.unwrapped)),s.push(Ze(l,t,n.linear)),o.push(Ze(l,t,n.residual))}return{ptfRaw:i,ptfUnwrapped:r,ptfLinear:s,ptfResidual:o}}function _s(n,t){const e=n.map((a,l)=>({dist:a,value:t[l]})).filter(a=>Number.isFinite(a.dist)&&Number.isFinite(a.value)).sort((a,l)=>a.dist-l.dist);if(e.length===0)return{dists:[],vals:[]};const i=Math.max(1,Math.min(16,Math.floor(e.length*.1)));let r=0,s=0;for(let a=0;a<i;a++)r+=e[a].value,s+=e[e.length-1-a].value;r/=i,s/=i;const o=r<=s?e:e.map(a=>({dist:-a.dist,value:a.value})).sort((a,l)=>a.dist-l.dist);return{dists:o.map(a=>a.dist),vals:o.map(a=>a.value)}}function mo(n,t,e="RGGB"){const i=e.toUpperCase(),r=t%2,s=n%2;return i==="RGGB"||i==="BGGR"?(r+s)%2!==0:i==="GBRG"||i==="GRBG"?(r+s)%2===0:(r+s)%2!==0}function go(n,t){return n+t&1?2:1}function Er(n,t){return(t&1)<<1|n&1}function yn(n,t,e){return n===void 0?0:typeof n=="number"?Number.isFinite(n)?n:0:Number.isFinite(n[Er(t,e)])?n[Er(t,e)]:0}function vt(n,t,e,i){return i!==void 0&&i!=="default"?go(n,t)===i:mo(n,t,e)}function ti(n){return n.length===0?0:n.reduce((t,e)=>t+e,0)/n.length}function yo(n,t,e){const i=(e%t+t)%t,r=Math.floor(i),s=(r+1)%t,o=i-r,a=r<n.length?n[r]:0,l=s<n.length?n[s]:0;return a*(1-o)+l*o}function Ms(n,t){const e=n.length;if(e===0)return 0;if(t<=0)return n[0];if(t>=e-1)return n[e-1];const i=Math.floor(t),r=Math.min(e-1,i+1),s=t-i;return n[i]*(1-s)+n[r]*s}function xo(n,t,e){const i=n.length,r=new Array(i).fill(0);for(let s=0;s<i;s++)r[s]=Ms(n,s-e+t);return r}function ei(n,t,e){const i=n.length;if(i===0)return{peakPos:0,peakIdx:0,peakVal:0};const r=Math.max(0,Math.floor(t-e)),s=Math.min(i-1,Math.ceil(t+e));let o=Math.max(0,Math.min(i-1,Math.round(t))),a=-1/0;for(let u=r;u<=s;u++){const d=Math.abs(n[u]);d>a&&(a=d,o=u)}Number.isFinite(a)||(a=Math.abs(n[o]??0));let l=o;if(o>0&&o<i-1){const u=n[o]>=0?1:-1,d=u*n[o-1],c=u*n[o],m=u*n[o+1],h=d-2*c+m;if(Number.isFinite(h)&&Math.abs(h)>1e-9){const f=.5*(d-m)/h;Number.isFinite(f)&&Math.abs(f)<=1&&(l=o+f)}}return{peakPos:l,peakIdx:o,peakVal:Math.abs(Ms(n,l))}}function bo(n,t,e,i,r,s,o,a,l){const u=Math.floor(i.x),d=Math.floor(i.y),c=Math.floor(i.w),m=Math.floor(i.h),h=[],f=(p,g)=>{if(p<0||g<0||p>=t||g>=e)return null;const x=r+p,y=s+g;return Math.max(0,n[g*t+p]-yn(l,x,y))};for(let p=0;p<m;p++){const g=[],x=d+p;for(let y=0;y<c;y++){const b=u+y,_=r+b,F=s+x;if(vt(_,F,o,a)){g.push(f(b,x)??0);continue}const v=[],w=f(b-1,x),P=f(b+1,x),T=f(b,x-1),C=f(b,x+1);if(w!==null&&vt(_-1,F,o,a)&&v.push(w),P!==null&&vt(_+1,F,o,a)&&v.push(P),T!==null&&vt(_,F-1,o,a)&&v.push(T),C!==null&&vt(_,F+1,o,a)&&v.push(C),v.length===0){const k=[],M=f(b-1,x-1),S=f(b+1,x-1),R=f(b-1,x+1),A=f(b+1,x+1);M!==null&&vt(_-1,F-1,o,a)&&k.push(M),S!==null&&vt(_+1,F-1,o,a)&&k.push(S),R!==null&&vt(_-1,F+1,o,a)&&k.push(R),A!==null&&vt(_+1,F+1,o,a)&&k.push(A),g.push(ti(k));continue}g.push(ti(v))}h.push(g)}return h}function ws(n,t,e,i,r,s,o,a,l){const u=Math.floor(i.x),d=Math.floor(i.y),c=Math.floor(i.w),m=Math.floor(i.h),h=[];for(let f=0;f<m;f++){const p=d+f;for(let g=0;g<c;g++){const x=u+g,y=r+x,b=s+p;vt(y,b,o,a)&&h.push({x:y,y:b,value:Math.max(0,n[p*t+x]-yn(l,y,b))})}}return h}function Ss(n,t,e,i,r){const s=Math.floor(i.x),o=Math.floor(i.y),a=Math.floor(i.w),l=Math.floor(i.h),u=(r==null?void 0:r.globalX)??0,d=(r==null?void 0:r.globalY)??0,c=!!(r!=null&&r.isThreePlane)&&n.length>=t*e*3,m=r==null?void 0:r.threePlaneChannel,h=[];for(let f=0;f<l;f++){const p=o+f,g=p*t;for(let x=0;x<a;x++){const y=s+x;let b=0;if(!c)b=Math.max(0,n[g+y]-yn(r==null?void 0:r.blackLevel,u+y,d+p));else{const _=(g+y)*3;if(m!==void 0)b=n[_+m];else{const F=n[_],v=n[_+1],w=n[_+2];b=.2126*F+.7152*v+.0722*w}}h.push({x:u+y,y:d+p,value:b})}}return h}function _o(n){var s;const t=n.length,e=((s=n[0])==null?void 0:s.length)??0;let i=0,r=0;for(let o=1;o<t-1;o++)for(let a=1;a<e-1;a++)i+=Math.abs(n[o][a+1]-n[o][a-1]),r+=Math.abs(n[o+1][a]-n[o-1][a]);return{gx:i,gy:r}}function vs(n,t,e,i,r,s,o){var m;const a=n.length,l=((m=n[0])==null?void 0:m.length)??0,u=(h,f,p,g,x,y,b)=>{const _=b?a:l,F=Math.max(0,f-3),v=Math.min(_,f+4);let w=0,P=0;for(let C=F;C<v;C++)w+=h[C],P+=C*h[C];if(w<=0)return null;const T=P/w;return b?{x:t+g*x,y:p+T*y,weight:w}:{x:p+T*y,y:e+g*x,weight:w}},d=(h,f,p,g,x,y,b)=>{const _=Math.max(3,Math.min(Math.max(3,Math.floor(p/3)),Math.max(4,Math.round(p*.12)))),F=h.map(S=>{let R=-1/0,A=-1;for(let L=0;L<S.length;L++)S[L]>R&&(R=S[L],A=L);return{peakValue:R,peakIndex:A}}),v=(f-1)*.5,w=F.map((S,R)=>({...S,index:R})).filter(S=>S.peakValue>1&&S.peakIndex>=0).sort((S,R)=>{const A=R.peakValue-S.peakValue;return Math.abs(A)>1e-6?A:Math.abs(S.index-v)-Math.abs(R.index-v)});if(w.length===0)return[];const P=w[0],T=new Array(f).fill(null),C=u(h[P.index],P.peakIndex,g,P.index,x,y,b);if(!C)return[];T[P.index]=C;const k=(S,R)=>{const A=h[S],L=F[S];if(!(L.peakValue>1)||L.peakIndex<0)return null;const X=Math.max(0,Math.floor(R-_)),D=Math.min(A.length,Math.ceil(R+_+1));let E=-1/0,U=-1;for(let V=X;V<D;V++)A[V]>E&&(E=A[V],U=V);if(U<0||!(E>1))return null;const O=Math.max(1e-6,L.peakValue),z=Math.abs(U-R)<=_,J=E>=O*.25;return!z||!J?null:u(A,U,g,S,x,y,b)};let M=C?b?(C.y-g)/y:(C.x-g)/y:P.peakIndex;for(let S=P.index+1;S<f;S++){const R=k(S,M);R&&(T[S]=R,M=b?(R.y-g)/y:(R.x-g)/y)}M=C?b?(C.y-g)/y:(C.x-g)/y:P.peakIndex;for(let S=P.index-1;S>=0;S--){const R=k(S,M);R&&(T[S]=R,M=b?(R.y-g)/y:(R.x-g)/y)}return T.filter(S=>!!S)};if(s){const h=n.map(f=>f.map((p,g)=>g===0?0:Math.abs(p-f[g-1])));return d(h,a,l,t,r,i,!1)}const c=[];for(let h=0;h<l;h++){const f=new Array(a).fill(0);for(let p=1;p<a;p++)f[p]=Math.abs(n[p][h]-n[p-1][h]);c.push(f)}return d(c,l,a,e,i,r,!0)}function Pe(n){if(n.length<2)return null;let t=0,e=0,i=0;for(const c of n)t+=c.weight,e+=c.x*c.weight,i+=c.y*c.weight;if(t<=0)return null;e/=t,i/=t;let r=0,s=0,o=0;for(const c of n){const m=c.x-e,h=c.y-i;r+=c.weight*m*m,s+=c.weight*h*h,o+=c.weight*m*h}r/=t,s/=t,o/=t;const a=.5*Math.atan2(2*o,r-s);let l=Math.cos(a),u=Math.sin(a);const d=Math.hypot(l,u);return!Number.isFinite(d)||d<=1e-9?null:(l/=d,u/=d,(l<0||Math.abs(l)<=1e-9&&u<0)&&(l=-l,u=-u),{pointX:e,pointY:i,dirX:l,dirY:u,orientation:Math.abs(l)>=Math.abs(u)?1:2})}function Mo(n,t){if(n.length!==4||t.length!==4||n.some(i=>i.length!==4))return null;const e=n.map((i,r)=>[...i,t[r]]);for(let i=0;i<4;i++){let r=i,s=Math.abs(e[i][i]);for(let a=i+1;a<4;a++){const l=Math.abs(e[a][i]);l>s&&(s=l,r=a)}if(!(s>1e-12))return null;if(r!==i){const a=e[i];e[i]=e[r],e[r]=a}const o=e[i][i];for(let a=i;a<=4;a++)e[i][a]/=o;for(let a=0;a<4;a++){if(a===i)continue;const l=e[a][i];if(!(Math.abs(l)<=1e-12))for(let u=i;u<=4;u++)e[a][u]-=l*e[i][u]}}return[e[0][4],e[1][4],e[2][4],e[3][4]]}function wo(n){if(n.length<4)return 0;const t=[...n].sort((m,h)=>m.x-h.x),e=t[0].x,r=t[t.length-1].x-e;if(!(r>1e-6))return 0;const s=16,o=[];for(let m=0;m<s;m++){const h=Math.max(0,Math.floor((m-1.5)*t.length/s)),f=Math.min(t.length-1,Math.floor((m+2.5)*t.length/s));if(f<h)continue;let p=0,g=0,x=0;for(let y=h;y<=f;y++)p+=t[y].x,g+=t[y].y,x++;x>0&&o.push({x:p/x,y:g/x})}if(o.length<4)return 0;const a=[.05952381,0,-.03571429,-.04761905,-.03571429,0,.05952381],l=new Array(o.length).fill(0),u=3;for(let m=0;m<o.length;m++){let h=0;for(let f=-u;f<=u;f++){const p=m+f;p<0||p>=o.length||(h+=a[f+u]*o[p].y)}l[m]=h}let d=0,c=1/0;for(let m=0;m<l.length-1;m++){const h=l[m],f=l[m+1];if(h===0){const b=Math.abs(o[m].x);b<c&&(c=b,d=o[m].x);continue}if(h*f>=0)continue;const p=f-h,g=Math.abs(p)>1e-12?-h/p:.5,x=o[m].x+(o[m+1].x-o[m].x)*g,y=Math.abs(x);y<c&&(c=y,d=x)}return!Number.isFinite(d)||d<e+.3*r||d>e+.7*r?0:d}function So(n){if(n.length<8)return null;const t=[...n].filter(f=>Number.isFinite(f.x)&&Number.isFinite(f.y)&&Number.isFinite(f.weight)).sort((f,p)=>f.x-p.x);if(t.length<8)return null;const i=[wo(t),0,.5*(t[Math.floor((t.length-1)*.5)].x+t[Math.ceil((t.length-1)*.5)].x)];let r=null;for(const f of i){if(!Number.isFinite(f))continue;let p=0,g=0;for(const x of t)x.x<=f?p++:g++;if(p>=4&&g>=4){r=f;break}}if(r===null)return null;const s=Array.from({length:4},()=>new Array(4).fill(0)),o=new Array(4).fill(0);for(const f of t){const p=f.x,g=f.y,x=Math.max(1e-6,f.weight),y=p<=r?[p*p,p,1,0]:[2*r*p-r*r,p,1,(p-r)*(p-r)];for(let b=0;b<4;b++){o[b]+=x*y[b]*g;for(let _=0;_<4;_++)s[b][_]+=x*y[b]*y[_]}}const a=Mo(s,o);if(!a)return null;const[l,u,d,c]=a,m=u+2*(l-c)*r,h=d+(c-l)*r*r;return[l,u,d,c,m,h].every(f=>Number.isFinite(f))?{splitX:r,left:[l,u,d],right:[c,m,h]}:null}function vo(n,t,e){const[i,r,s]=e;if(Math.abs(i)<=1e-12){const b=1+r*r;return b>1e-12?[(n-r*(s-t))/b]:[n]}const o=2*i*i,a=3*i*r,l=1+2*i*s-2*i*t+r*r,u=r*s-t*r-n;if(Math.abs(o)<=1e-12)return[n];const d=a/o,c=l/o,m=u/o,h=(d*d-3*c)/9,f=(2*d*d*d-9*d*c+27*m)/54,p=f*f-h*h*h;if(p<0&&h>0){const b=Math.acos(Math.max(-1,Math.min(1,f/Math.sqrt(h*h*h)))),_=-2*Math.sqrt(h);return[_*Math.cos(b/3)-d/3,_*Math.cos((b+2*Math.PI)/3)-d/3,_*Math.cos((b-2*Math.PI)/3)-d/3]}const g=Math.sqrt(Math.max(0,p)),x=-Math.sign(f||1)*Math.cbrt(Math.abs(f)+g),y=Math.abs(x)<=1e-12?0:h/x;return[x+y-d/3]}function Ps(n,t){if(n.length<8)return null;const e=-t.dirY,i=t.dirX,r=n.map(o=>({x:(o.x-t.pointX)*t.dirX+(o.y-t.pointY)*t.dirY,y:(o.x-t.pointX)*e+(o.y-t.pointY)*i,weight:o.weight})),s=So(r);return s?{...t,normalX:e,normalY:i,splitX:s.splitX,left:s.left,right:s.right}:null}function Po(n,t){const e=n.x-t.pointX,i=n.y-t.pointY,r=e*t.dirX+i*t.dirY,s=e*t.normalX+i*t.normalY,o=r<t.splitX?t.left:t.right,a=vo(r,s,o);let l=s,u=Number.POSITIVE_INFINITY;for(const d of a){if(!Number.isFinite(d))continue;const c=o[0]*d*d+o[1]*d+o[2],m=r-d,h=s-c,f=Math.hypot(m,h);Number.isFinite(f)&&f<u&&(u=f,l=(h>=0?1:-1)*f)}return Number.isFinite(u)?l:s}function Cs(n,t,e,i,r,s,o,a=Lt){if(!t||t.length<8||n.length===0)return null;const l=t.filter(S=>Number.isFinite(S.x)&&Number.isFinite(S.y)).map(S=>({x:S.x,y:S.y,weight:1}));if(l.length<8)return null;const u=Pe(l);if(!u)return null;const d=e.p2.x-e.p1.x,c=e.p2.y-e.p1.y,m=Math.hypot(d,c);if(!Number.isFinite(m)||m<=1e-6)return null;let h=u.dirX,f=u.dirY;h*d+f*c<0&&(h=-h,f=-f);const p={...u,dirX:h,dirY:f},g=Ps(l,p),x=d/m,y=c/m,b=-y,_=x,F=(e.p1.x+e.p2.x)*.5,v=(e.p1.y+e.p2.y)*.5,w=-p.dirY,P=p.dirX,T=Math.abs(x)>=Math.abs(y)?1:2,C=[],k=[];for(const S of n){const R=S.x-F,A=S.y-v,L=R*x+A*y;if(Math.abs(L)>i)continue;const X=R*b+A*_;if(Math.abs(X)>r)continue;const D=g?Po(S,g):(S.x-p.pointX)*w+(S.y-p.pointY)*P;Number.isFinite(D)&&(C.push(D),k.push(S.value))}if(C.length<8)return null;const M=o?s!=null&&s.forceLegacyModel?Fn(C,k,T,a,r*2):xn(C,k,T,a):ui(C,k,Math.max(2,r*2),s==null?void 0:s.manualBinSize,T,s==null?void 0:s.preferAutoPerEdgeBin);return M?(M.quadraticProjectionUsed=!!g,M):null}function Ur(n){if(n.length<2)return null;const t=n.filter(e=>Number.isFinite(e.x)&&Number.isFinite(e.y)).map(e=>({x:e.x,y:e.y,weight:1}));return t.length<2?null:vn(t,Pe(t))}function vn(n,t){if(!t||n.length<2)return null;let e=1/0,i=-1/0;for(const s of n){const o=(s.x-t.pointX)*t.dirX+(s.y-t.pointY)*t.dirY;e=Math.min(e,o),i=Math.max(i,o)}if(!Number.isFinite(e)||!Number.isFinite(i))return null;const r=Math.max(.5,(i-e)*.03);return{p1:{x:t.pointX+t.dirX*(e-r),y:t.pointY+t.dirY*(e-r)},p2:{x:t.pointX+t.dirX*(i+r),y:t.pointY+t.dirY*(i+r)}}}function Co(n,t,e,i,r,s,o,a){return[vn(n,t),vn(e,i),vn(r,s),vn(o,a)]}function ko(n,t,e){if(!n||n.length<8)return;const i=n.filter(y=>Number.isFinite(y.x)&&Number.isFinite(y.y)).map(y=>({x:y.x,y:y.y,weight:1}));if(i.length<8)return;const r=Pe(i);if(!r)return;const s=t.p2.x-t.p1.x,o=t.p2.y-t.p1.y,a=Math.hypot(s,o);if(!Number.isFinite(a)||a<=1e-6)return;let l=r.dirX,u=r.dirY;l*s+u*o<0&&(l=-l,u=-u);const d=Ps(i,{...r,dirX:l,dirY:u});if(!d)return;const c=n.map(y=>(y.x-d.pointX)*d.dirX+(y.y-d.pointY)*d.dirY).filter(y=>Number.isFinite(y)),m=(t.p1.x-d.pointX)*d.dirX+(t.p1.y-d.pointY)*d.dirY,h=(t.p2.x-d.pointX)*d.dirX+(t.p2.y-d.pointY)*d.dirY;if(Number.isFinite(m)&&c.push(m),Number.isFinite(h)&&c.push(h),c.length<2)return;const f=Math.min(...c),p=Math.max(...c);if(!Number.isFinite(f)||!Number.isFinite(p)||p-f<=1e-6)return;const g=Math.max(21,e),x=[];for(let y=0;y<g;y++){const b=g===1?.5:y/(g-1),_=f+(p-f)*b,F=_<d.splitX?d.left:d.right,v=F[0]*_*_+F[1]*_+F[2];x.push({x:d.pointX+_*d.dirX+v*d.normalX,y:d.pointY+_*d.dirY+v*d.normalY})}return x}function ui(n,t,e,i,r,s=!1,o=!1,a=!1){if(n.length===0||t.length!==n.length)return null;const l=_s(n,t),u=l.dists,d=l.vals;if(u.length===0)return null;const c=()=>{const w=e/2;let P=0;for(const C of n)Math.abs(C)<=w&&P++;if(P<=0)return .125;const T=40*w/P;return Math.max(.01,Math.min(.125,T))},m=(w,P,T,C,k)=>{if(!(C>0)||!(k>0)||!(T>P))return!1;const M=Math.floor((T-P)/C);if(M<2)return!1;const S=Math.max(P,-k),R=Math.min(T,k);if(!(R>S))return!1;const A=Math.max(0,Math.floor((S-P)/C)),L=Math.min(M,Math.ceil((R-P)/C));if(L<=A)return!1;const X=new Array(L-A).fill(0),D=P+A*C,E=P+L*C;for(let U=0;U<w.length;U++){const O=w[U];if(O<D)continue;if(O>=E)break;const z=Math.floor((O-P)/C);z>=A&&z<L&&X[z-A]++}return X.every(U=>U>0)},h=()=>{const w=u[0],P=u[u.length-1],T=Math.max(0,e*.25),C=.125,k=.5,M=.001,S=Math.round((k-C)/M);for(let R=0;R<=S;R++){const A=Number((C+R*M).toFixed(3));if(m(u,w,P,A,T))return A}return k};let f=.125;i&&i>0?f=Math.max(.01,Math.min(.5,i)):s?f=h():f=c();const p=u[0],g=u[u.length-1],x=Math.floor((g-p)/f);if(x<2)return null;const y=()=>{const w=new Array(x).fill(0),P=new Array(x).fill(0);for(let C=0;C<u.length;C++){const k=(u[C]-p)/f;if(Number.isFinite(k))if(o){const M=Math.floor(k),S=k-M,R=1-S,A=S;M>=0&&M<x&&(w[M]+=d[C]*R,P[M]+=R);const L=M+1;L>=0&&L<x&&(w[L]+=d[C]*A,P[L]+=A)}else{const M=Math.floor(k);M>=0&&M<x&&(w[M]+=d[C],P[M]++)}}let T=d[0];for(let C=0;C<x;C++)P[C]>0?(w[C]/=P[C],T=w[C]):w[C]=T;return w},b=a?null:oo(u,d,f),_=(b==null?void 0:b.esf)??y(),F=(b==null?void 0:b.lsfFull)??no(_,f),v=Math.max(0,Math.min(x-1,-p/f-.5));return{esf:_,lsfFull:F,binSize:f,orientation:r,zeroIndex:v,shortSidePx:e,fallbackUsed:a||!b}}function xn(n,t,e,i=Lt){if(n.length===0||t.length!==n.length)return null;const r=_s(n,t),s=r.dists.map((a,l)=>({dist:a,value:r.vals[l]})).filter(a=>Math.abs(a.dist)<i);if(s.length<8)return null;const o=ys(s.map(a=>a.dist),s.map(a=>a.value),i);return o?{esf:o.esf,lsfFull:o.lsfFull,binSize:ci,orientation:e,zeroIndex:ae/2,shortSidePx:i*2,fallbackUsed:!1,mtfmapperLike:!0,mtfmapperOrderedDists:s.map(a=>a.dist),mtfmapperOrderedVals:s.map(a=>a.value),mtfmapperEffectiveMaxDot:i}:null}function Fn(n,t,e,i=Lt,r=i*2){if(n.length===0||t.length!==n.length)return null;const s=[],o=[];for(let a=0;a<n.length;a++){const l=n[a],u=t[a];!Number.isFinite(l)||!Number.isFinite(u)||Math.abs(l)>=i||(s.push(l),o.push(u))}return s.length<8?null:ui(s,o,Math.max(2,r),void 0,e,!0,!0,!0)}function ve(n,t,e,i,r,s,o,a,l=Lt){if(s<=0||o<=0)return null;const d=!(!!(a!=null&&a.isThreePlane)&&n.length>=t*e*3)&&((a==null?void 0:a.greenOnly)??!1),c=(a==null?void 0:a.bayerPattern)||"RGGB",m=r.p2.x-r.p1.x,h=r.p2.y-r.p1.y,f=Math.hypot(m,h);if(!Number.isFinite(f)||f<=1e-6)return null;const p=m/f,g=h/f,x=-g,y=p,b=(r.p1.x+r.p2.x)*.5,_=(r.p1.y+r.p2.y)*.5,F=Math.abs(p)>=Math.abs(g)?1:2,v=d?ws(n,t,e,i,0,0,c,a==null?void 0:a.greenPhase,a==null?void 0:a.blackLevel):Ss(n,t,e,i,{...a,globalX:0,globalY:0});if(v.length===0)return null;if(!(a!=null&&a.disableQuadraticProjection)){const T=Cs(v,a==null?void 0:a.quadraticFitPoints,r,s,o,a,!0,l);if(T)return T}const w=[],P=[];for(const T of v){const C=T.x-b,k=T.y-_,M=C*p+k*g;if(Math.abs(M)>s)continue;const S=C*x+k*y;Math.abs(S)>o||(w.push(S),P.push(T.value))}return w.length<8?null:a!=null&&a.forceLegacyModel?Fn(w,P,F,l):xn(w,P,F,l)}function Fo(n,t,e=0){const i=[...n.lsfFull];if(i.length<3)return!1;const r=Math.max(n.binSize,1e-6),s=Number.isFinite(n.zeroIndex)?n.zeroIndex:i.length/2,o=Math.max(1,Math.round((n.shortSidePx??0)*.5/r));let{peakPos:a,peakIdx:l,peakVal:u}=ei(i,s,o);const d=u*.2;let c=0,m=i.length-1;for(let p=l;p>=0;p--)if(i[p]<d){c=p;break}for(let p=l;p<i.length;p++)if(i[p]<d){m=p;break}const h=m-c;if(t&&h>0){const p=h*4,g=[],x=[];if(e>0){const y=Math.max(0,l-p-e),b=Math.max(0,l-p);for(let v=y;v<b;v++)g.push(v),x.push(i[v]);const _=Math.min(i.length,l+p),F=Math.min(i.length,l+p+e);for(let v=_;v<F;v++)g.push(v),x.push(i[v])}else{for(let y=0;y<Math.max(0,l-p);y++)g.push(y),x.push(i[y]);for(let y=Math.min(i.length,l+p);y<i.length;y++)g.push(y),x.push(i[y])}if(g.length>2){const{slope:y,intercept:b}=Zi(g,x);for(let _=0;_<i.length;_++)i[_]=i[_]-(y*_+b);({peakPos:a}=ei(i,s,o))}}return Math.abs(a-s)*r<=Math.max(1e-6,(n.shortSidePx??0)/6)}function Ao(n){const t=n.length;if(t<3)return!1;let e=0,i=-1/0;for(let o=0;o<t;o++){const a=Math.abs(n[o]);a>i&&(i=a,e=o)}const r=t/3,s=2*t/3;return e>=r&&e<=s}function Mt(n,t,e){const i=Math.max(0,Math.floor(n.x)),r=Math.max(0,Math.floor(n.y)),s=Math.min(t,Math.ceil(n.x+n.w)),o=Math.min(e,Math.ceil(n.y+n.h)),a=s-i,l=o-r;return a<2||l<2?null:{x:i,y:r,w:a,h:l}}function er(n,t,e,i){const r=[],s=n.x,o=n.y,a=n.x+n.w,l=n.y+n.h,u=n.x+n.w*.5,d=n.y+n.h*.5,c=[{x:s,y:o},{x:a,y:o},{x:a,y:l},{x:s,y:l},{x:u,y:o},{x:a,y:d},{x:u,y:l},{x:s,y:d},{x:u,y:d}];for(const m of c){const h=ke(m,t);Number.isFinite(h.x)&&Number.isFinite(h.y)&&r.push(h)}return r.length===0?null:Mt(Vt(r,2),e,i)}function Vt(n,t=0){let e=1/0,i=1/0,r=-1/0,s=-1/0;for(const o of n)e=Math.min(e,o.x),i=Math.min(i,o.y),r=Math.max(r,o.x),s=Math.max(s,o.y);return{x:e-t,y:i-t,w:r-e+t*2,h:s-i+t*2}}function Dr(n,t){let e=Math.atan2(t,n)*180/Math.PI;return e<0&&(e+=180),e}function be(n,t){const e=n.p2.x-n.p1.x,i=n.p2.y-n.p1.y,r=Math.hypot(e,i);if(!Number.isFinite(r)||r<=1e-6)return null;const s=-i/r,o=e/r;return[{x:n.p1.x+s*t,y:n.p1.y+o*t},{x:n.p2.x+s*t,y:n.p2.y+o*t},{x:n.p2.x-s*t,y:n.p2.y-o*t},{x:n.p1.x-s*t,y:n.p1.y-o*t}]}function To(n,t,e,i,r,s,o,a){if(s<=0||o<=0)return null;const u=!(!!(a!=null&&a.isThreePlane)&&n.length>=t*e*3)&&((a==null?void 0:a.greenOnly)??!1),d=(a==null?void 0:a.bayerPattern)||"RGGB",c=r.p2.x-r.p1.x,m=r.p2.y-r.p1.y,h=Math.hypot(c,m);if(!Number.isFinite(h)||h<=1e-6)return null;const f=c/h,p=m/h,g=-p,x=f,y=(r.p1.x+r.p2.x)*.5,b=(r.p1.y+r.p2.y)*.5,_=Math.abs(f)>=Math.abs(p)?1:2,F=u?ws(n,t,e,i,0,0,d,a==null?void 0:a.greenPhase,a==null?void 0:a.blackLevel):Ss(n,t,e,i,{...a,globalX:0,globalY:0});if(F.length===0)return null;if(!(a!=null&&a.disableQuadraticProjection)){const P=Cs(F,a==null?void 0:a.quadraticFitPoints,r,s,o,a,!1);if(P)return P}const v=[],w=[];for(const P of F){const T=P.x-y,C=P.y-b,k=T*f+C*p;if(Math.abs(k)>s)continue;const M=T*g+C*x;Math.abs(M)>o||(v.push(M),w.push(P.value))}return v.length<8?null:ui(v,w,o*2,a==null?void 0:a.manualBinSize,_,a==null?void 0:a.preferAutoPerEdgeBin)}function Io(n,t,e){const i=[...n],r=new Array(n.length).fill(0),s=[0,0,0];let o=-1,a=1,l=-1;for(let c=1;c<n.length-1;c++){let m=0;if(n[c]>1e-4){m=Math.atan2(e[c]*o,t[c]*o);let h=0;for(let f=-5;f<=5;f++)Math.abs(m+f*2*Math.PI-s[1])<Math.abs(m+h*2*Math.PI-s[1])&&(h=f);m+=h*2*Math.PI}c>3&&Math.abs(m-s[0])>Math.PI/2&&l<c-1&&n[c]<.5&&(a*=-1,l=c),i[c]*=a,o*=-1,s[0]=s[1],s[1]=m,s[2]=m}const u=[-.086,.343,.486,.343,-.086];for(let c=0;c<n.length-3;c++){let m=0;for(let h=-2;h<=2;h++)m+=i[Math.abs(c+h)]*u[h+2];r[c]=m}for(let c=0;c<n.length-3;c++)i[c]=r[c];const d=7;for(let c=0;c<3;c++){r.fill(0);for(let h=0;h<n.length-d;h++)if(h<d)r[h]=i[h];else{const f=Math.min(5,Math.floor((h-5)/3)),p=io[f];let g=0;for(let x=-d;x<=d;x++)g+=i[h+x]*p[x+d];r[h]=g}for(let h=n.length-d-2;h<n.length;h++)r[h]=i[h];const m=Math.abs(r[0])>1e-9?r[0]:1;for(let h=0;h<n.length;h++)i[h]=r[h]/m}for(let c=0;c<n.length;c++)i[c]=Math.abs(i[c]);return i}function Ro(n,t){const e=[[0,0,0],[0,0,0],[0,0,0]],i=[0,0,0];for(let o=0;o<n.length;o++){const a=n[o],l=-t+o,u=[1,a,a*a];for(let d=0;d<3;d++){i[d]+=u[d]*l;for(let c=0;c<3;c++)e[d][c]+=u[d]*u[c]}}const r=Gs(e);if(!r)return null;const s=Xs(r,i);return[s[0],s[1],s[2]]}function No(n,t){let e=0,i=1,r=0,s=!1,o=0;const a=Math.min(n.length,ae/16*2);for(let l=0;l<a&&!s;l++){const u=n[l];if(i>.5&&u<=.5){const d=-(u-i)*ae;Math.abs(d)>1e-9&&(r=-((.5-i-d*e)/d),o=l,s=!0)}i=u,e=l/ae}if(!s)return null;if(o>=5&&o<a-10){const l=Math.min(Math.max(2,o-9),9),u=Ro(n.slice(o-l,o+l+1),l);if(u){const c=(u[0]+.5*u[1]+.25*u[2]+o)/ae;if(o>9)r=c;else{const h=(o-5)/ae/8;r=(1-h)*r+h*c}}}return r*Et*t}function Lo(n,t){if(n.length===0)return null;const e=ae,i=ae/16*4,r=new kn(e),s=1,o=ao(),a=new Float32Array(501);for(let M=0;M<=500;M++)a[M]=M/500*s*2;const l=new Array(i).fill(0).map((M,S)=>S/e*s*Et),u=new Float32Array(i).fill(0),d=new Float32Array(i).fill(0);let c=0,m=[],h=[],f=[],p=[],g=[],x=[],y=[],b=null,_=0;for(const M of n){const S=M.mtfmapperOrderedDists&&M.mtfmapperOrderedVals&&M.mtfmapperOrderedDists.length===M.mtfmapperOrderedVals.length?ys(M.mtfmapperOrderedDists,M.mtfmapperOrderedVals,M.mtfmapperEffectiveMaxDot??Lt):null,R=(S==null?void 0:S.lsfFull)??M.lsfFull,A=(S==null?void 0:S.esf)??M.esf;if(R.length<e)continue;const L=new Float32Array(e);for(let E=0;E<e;E++)L[E]=R[E]??0;r.transform(L);const X=Math.max(1e-9,Math.abs(r._real[0])),D=new Array(i).fill(0);for(let E=1;E<i;E++)D[E]=Math.atan2(r._imag[E],r._real[E]);for(let E=0;E<i;E++)u[E]+=r._real[E]/X,d[E]+=r._imag[E]/X;if(c++,_+=M.shortSidePx*.5,m.length===0){m=[...R],h=[...A];const E=new Array(i).fill(0);E[0]=1;for(let V=1;V<i;V++)E[V]=Math.hypot(r._real[V]/X,r._imag[V]/X);const U=l.map(V=>V),O=(Number.isFinite(M.zeroIndex)?M.zeroIndex:0)*(M.binSize??ci),z=xs(D,U,E,Number.POSITIVE_INFINITY,O),J=bs(z,l,a);p=J.ptfRaw,g=J.ptfUnwrapped,x=J.ptfLinear,y=J.ptfResidual,f=J.ptfResidual,b=z.fit}}if(c===0)return null;const F=new Float32Array(i),v=new Float32Array(i),w=new Array(i).fill(0);w[0]=1;for(let M=0;M<i;M++)F[M]=u[M]/c,v[M]=d[M]/c,M>0&&(w[M]=Math.hypot(F[M],v[M]));const P=Io(w,F,v),T=new Array(i).fill(0);for(let M=0;M<i;M++)T[M]=P[M]/o[M];const C=Array.from(a,M=>Ze(M,l,T)),k=No(T,s);return{esf:h,lsf:[],lsfCropped:m,mtf:C,ptf:f,ptfRaw:p,ptfUnwrapped:g,ptfLinear:x,ptfResidual:y,ptfPhaseFit:b,freqs:Array.from(a),mtf50:k,calcRadius:_/c}}function Eo(n,t,e,i=!1,r=0,s=!1){if(n.length===0)return null;if(n.every(w=>w.mtfmapperLike))return Lo(n);const o=4096,a=new kn(o),l=1,u=new Float32Array(501);for(let w=0;w<=500;w++)u[w]=w/500*l*2;const d=new Float32Array(501).fill(0);let c=0,m=[],h=[],f=0,p=[],g=[],x=[],y=[],b=[],_=null;for(const w of n){let P=[...w.lsfFull];const T=w.binSize,C=Number.isFinite(w.zeroIndex)?w.zeroIndex:P.length/2,k=Math.max(1,Math.round((w.shortSidePx??0)*.5/Math.max(T,1e-6)));let{peakPos:M,peakIdx:S,peakVal:R}=ei(P,C,k);const A=R*.2;let L=0,X=P.length-1;for(let I=S;I>=0;I--)if(P[I]<A){L=I;break}for(let I=S;I<P.length;I++)if(P[I]<A){X=I;break}const D=X-L;let E=!1;if(i&&D>0){const I=D*4,H=[],G=[];if(r>0){const Y=Math.max(0,S-I-r),et=Math.max(0,S-I);for(let ct=Y;ct<et;ct++)H.push(ct),G.push(P[ct]);const nt=Math.min(P.length,S+I),it=Math.min(P.length,S+I+r);for(let ct=nt;ct<it;ct++)H.push(ct),G.push(P[ct])}else{for(let Y=0;Y<Math.max(0,S-I);Y++)H.push(Y),G.push(P[Y]);for(let Y=Math.min(P.length,S+I);Y<P.length;Y++)H.push(Y),G.push(P[Y])}if(H.length>2){const{slope:Y,intercept:et}=Zi(H,G);for(let nt=0;nt<P.length;nt++)P[nt]=P[nt]-(Y*nt+et);({peakPos:M,peakIdx:S,peakVal:R}=ei(P,C,k)),E=!0}}let U=0,O=0;if(t>0)O=t,U=Math.round(t/T);else{const I=R*.2;let H=0,G=P.length-1;for(let it=S;it>=0;it--)if(P[it]<I){H=it;break}for(let it=S;it<P.length;it++)if(P[it]<I){G=it;break}const et=(G-H)*T;let nt=Math.max(2,et*8);O=nt,U=Math.round(nt/T)}f+=O;const z=Math.max(0,Math.floor(C-U)),J=Math.min(P.length,Math.ceil(C+U)),V=P.slice(z,J);if(V.length<8)continue;const Q=new Float32Array(o).fill(0),W=new Array(V.length).fill(0);for(let I=0;I<V.length;I++){let H=1;s&&(H=.5*(1-Math.cos(2*Math.PI*I/(V.length-1)))),W[I]=V[I]*H}const tt=Math.max(0,Math.min(V.length-1,M-z));for(let I=0;I<o;I++)Q[I]=yo(W,o,I+tt);a.transform(Q);const at=[],B=[],Z=[];for(let I=0;I<=o/2;I++){const H=a._real[I],G=a._imag[I],Y=Math.sqrt(H*H+G*G);at.push(Y),B.push(I/(o*T)*l),Z.push(Math.atan2(G,H))}const q=at[0];if(q>0){for(let I=0;I<=500;I++){const H=u[I],Y=so(H,T);d[I]+=Ze(H,B,at)/q/Y}if(c++,m.length===0){m=xo(V,tt,(V.length-1)/2),h=E?Uo(P):w.esf;const I=at.map(et=>et/q),H=B.map(et=>et),G=xs(Z,H,I,Number.POSITIVE_INFINITY,0),Y=bs(G,B,u);g=Y.ptfRaw,x=Y.ptfUnwrapped,y=Y.ptfLinear,b=Y.ptfResidual,p=Y.ptfResidual,_=G.fit}}}if(c===0)return null;const F=Array.from(d).map(w=>w/c);let v=null;for(let w=0;w<F.length-1;w++)if(F[w]>=.5&&F[w+1]<.5){v=u[w]+(.5-F[w])*(u[w+1]-u[w])/(F[w+1]-F[w]);break}return{esf:h,lsf:[],lsfCropped:m,mtf:F,ptf:p,ptfRaw:g,ptfUnwrapped:x,ptfLinear:y,ptfResidual:b,ptfPhaseFit:_,freqs:Array.from(u),mtf50:v,calcRadius:f/c}}function Uo(n){const t=new Array(n.length).fill(0);let e=0;for(let i=0;i<n.length;i++)e+=n[i],t[i]=e;return t}function Ze(n,t,e){if(n<=t[0])return e[0];if(n>=t[t.length-1])return e[e.length-1];let i=0;for(;n>t[i+1];)i++;const r=(n-t[i])/(t[i+1]-t[i]);return e[i]+r*(e[i+1]-e[i])}function nr(n){return{...zi,...n,gradientPercentiles:n!=null&&n.gradientPercentiles&&n.gradientPercentiles.length>0?n.gradientPercentiles:zi.gradientPercentiles}}function Do(n){return!n||n.length===0?void 0:[Number.isFinite(n[0])?n[0]:0,Number.isFinite(n[1])?n[1]:Number.isFinite(n[0])?n[0]:0,Number.isFinite(n[2])?n[2]:Number.isFinite(n[0])?n[0]:0,Number.isFinite(n[3])?n[3]:Number.isFinite(n[0])?n[0]:0]}function Bo(n,t){const e=n.width,i=n.height,r=n.data,s=n.bayerPattern||"RGGB",o=Do(n.blackLevels),a=new Float32Array(e*i),l=(_,F)=>_<0||F<0||_>=e||F>=i?null:Math.max(0,r[F*e+_]-yn(o,_,F));let u=1/0,d=-1/0;for(let _=0;_<i;_++){const F=_*e;for(let v=0;v<e;v++){const w=F+v;let P=0;if(vt(v,_,s,t))P=l(v,_)??0;else{const T=[],C=l(v-1,_),k=l(v+1,_),M=l(v,_-1),S=l(v,_+1);if(C!==null&&vt(v-1,_,s,t)&&T.push(C),k!==null&&vt(v+1,_,s,t)&&T.push(k),M!==null&&vt(v,_-1,s,t)&&T.push(M),S!==null&&vt(v,_+1,s,t)&&T.push(S),T.length>0)P=ti(T);else{const R=[],A=l(v-1,_-1),L=l(v+1,_-1),X=l(v-1,_+1),D=l(v+1,_+1);A!==null&&vt(v-1,_-1,s,t)&&R.push(A),L!==null&&vt(v+1,_-1,s,t)&&R.push(L),X!==null&&vt(v-1,_+1,s,t)&&R.push(X),D!==null&&vt(v+1,_+1,s,t)&&R.push(D),P=ti(R)}}a[w]=P,P<u&&(u=P),P>d&&(d=P)}}if(!Number.isFinite(u)||!Number.isFinite(d)||d<=u+1e-9)return new Uint8Array(e*i);const c=1024,m=new Uint32Array(c),h=d-u;for(let _=0;_<a.length;_++){const F=Math.max(0,Math.min(1,(a[_]-u)/h)),v=Math.min(c-1,Math.max(0,Math.floor(F*(c-1))));m[v]++}const f=a.length,p=_=>{const F=f*_;let v=0;for(let w=0;w<c;w++)if(v+=m[w],v>=F)return u+w/Math.max(1,c-1)*h;return d},g=p(.01),x=p(.99),y=Math.max(1e-9,x-g),b=new Uint8Array(e*i);for(let _=0;_<a.length;_++){const F=Math.max(0,Math.min(1,(a[_]-g)/y));b[_]=Math.round(F*255)}return b}function Oo(n,t,e){const i=new Float32Array(n.length),r=new Float32Array(n.length),s=new Float32Array(n.length);for(let o=1;o<e-1;o++)for(let a=1;a<t-1;a++){const l=o*t+a,u=n[(o-1)*t+(a-1)],d=n[(o-1)*t+a],c=n[(o-1)*t+(a+1)],m=n[o*t+(a-1)],h=n[o*t+(a+1)],f=n[(o+1)*t+(a-1)],p=n[(o+1)*t+a],g=n[(o+1)*t+(a+1)],x=-u-2*m-f+(c+2*h+g),y=-u-2*d-c+(f+2*p+g);i[l]=x,r[l]=y,s[l]=Math.hypot(x,y)}return{gx:i,gy:r,magnitude:s}}function Vo(n,t){let e=0,i=0;for(let l=0;l<n.length;l++){const u=n[l];!Number.isFinite(u)||u<=1e-6||(e=Math.max(e,u),i++)}if(i===0||e<=1e-6)return[];const r=1024,s=new Uint32Array(r);for(let l=0;l<n.length;l++){const u=n[l];if(!Number.isFinite(u)||u<=1e-6)continue;const d=Math.max(0,Math.min(1,u/e)),c=Math.min(r-1,Math.floor(d*(r-1)));s[c]++}const o=t&&t.length>0?t:zi.gradientPercentiles,a=[];for(const l of o){const u=i*l;let d=0;for(let c=0;c<r;c++)if(d+=s[c],d>=u){a.push(c/Math.max(1,r-1)*e);break}}return Array.from(new Set(a.filter(l=>l>0))).sort((l,u)=>u-l)}function Go(n,t){const e=new Uint8Array(n.length);for(let i=0;i<n.length;i++)e[i]=n[i]>=t?1:0;return e}const Xo=256*256;function zo(n,t,e){if(n.length>=Xo){const s=to.compute(n,t,e);if(s)return{gray:s.blurredGray,gradient:{gx:s.gx,gy:s.gy,magnitude:s.magnitude},backend:"webgl"}}const i=Zo(n,t,e),r=Oo(i,t,e);return{gray:i,gradient:r,backend:"cpu"}}function Yo(n,t,e,i){let r=n;for(let s=0;s<i;s++){const o=new Uint8Array(n.length);for(let a=0;a<e;a++)for(let l=0;l<t;l++){let u=0;for(let d=-1;d<=1&&!u;d++){const c=a+d;if(!(c<0||c>=e))for(let m=-1;m<=1;m++){const h=l+m;if(!(h<0||h>=t)&&r[c*t+h]){u=1;break}}}o[a*t+l]=u}r=o}return r}function Wo(n,t,e){const i=new Int32Array(n.length),r=[];let s=1;for(let o=0;o<n.length;o++){if(!n[o]||i[o]!==0)continue;const a=[o];i[o]=s;let l=0,u=t,d=e,c=0,m=0,h=0,f=!1;for(;l<a.length;){const p=a[l++],g=p%t,x=Math.floor(p/t);h++,u=Math.min(u,g),d=Math.min(d,x),c=Math.max(c,g),m=Math.max(m,x),(g===0||x===0||g===t-1||x===e-1)&&(f=!0);for(let y=-1;y<=1;y++)for(let b=-1;b<=1;b++){if(b===0&&y===0)continue;const _=g+b,F=x+y;if(_<0||F<0||_>=t||F>=e)continue;const v=F*t+_;!n[v]||i[v]!==0||(i[v]=s,a.push(v))}}r.push({label:s,x:u,y:d,w:c-u+1,h:m-d+1,area:h,touchesBorder:f}),s++}return{labels:i,components:r}}function ks(n,t){const e=Math.hypot(n,t);if(!Number.isFinite(e)||e<=1e-9)return null;let i=n/e,r=t/e;return(i<0||Math.abs(i)<=1e-9&&r<0)&&(i=-i,r=-r),{x:i,y:r}}function Te(n,t){if(n.length===0)return 0;const e=[...n].sort((o,a)=>o.value-a.value),i=e.reduce((o,a)=>o+Math.max(0,a.weight),0);if(i<=0)return e[Math.floor((e.length-1)*t)].value;const r=Math.max(0,Math.min(1,t))*i;let s=0;for(const o of e)if(s+=Math.max(0,o.weight),s>=r)return o.value;return e[e.length-1].value}function Br(n){const t=n.filter(i=>Number.isFinite(i)).sort((i,r)=>i-r);if(t.length===0)return 0;const e=i=>{if(i.length===1)return i[0];if(i.length===2)return(i[0]+i[1])*.5;const r=Math.ceil(i.length*.5);let s=0,o=1/0;for(let a=0;a+r-1<i.length;a++){const l=i[a+r-1]-i[a];l<o&&(o=l,s=a)}return e(i.slice(s,s+r))};return e(t)}function jo(n,t,e,i,r,s,o){const a=[];for(let l=r.y;l<r.y+r.h;l++)for(let u=r.x;u<r.x+r.w;u++){const d=l*s+u;if(n[d]!==t||!e[d])continue;const c=i.magnitude[d];!Number.isFinite(c)||c<=1e-6||a.push({x:u,y:l,weight:c,gx:i.gx[d],gy:i.gy[d]})}return a}function Ho(n){let t=0,e=0,i=0,r=0,s=0;for(const l of n){t+=l.weight,e+=l.x*l.weight,i+=l.y*l.weight;const u=Math.hypot(l.gx,l.gy);if(!Number.isFinite(u)||u<=1e-6)continue;const d=-l.gy/u,c=l.gx/u;r+=l.weight*(d*d-c*c),s+=l.weight*(2*d*c)}if(t<=0)return null;e/=t,i/=t;const o=.5*Math.atan2(s,r),a=ks(Math.cos(o),Math.sin(o));return a?{centerX:e,centerY:i,dirX:a.x,dirY:a.y,orthoX:-a.y,orthoY:a.x}:null}function Yn(n,t){let e=0,i=0;const r=-t.dirY,s=t.dirX;for(const o of n){const a=(o.x-t.pointX)*r+(o.y-t.pointY)*s;i+=o.weight*a*a,e+=o.weight}return e<=0?1/0:Math.sqrt(i/e)}function Or(n,t,e,i,r){const s=Math.max(0,Math.min(t-1,i)),o=Math.max(0,Math.min(e-1,r)),a=Math.floor(s),l=Math.floor(o),u=Math.min(t-1,a+1),d=Math.min(e-1,l+1),c=s-a,m=o-l,h=n[l*t+a],f=n[l*t+u],p=n[d*t+a],g=n[d*t+u],x=h+(f-h)*c,y=p+(g-p)*c;return x+(y-x)*m}function Qo(n,t,e,i,r,s,o,a){const l=Or(n,t,e,i-s*a,r-o*a);return Or(n,t,e,i+s*a,r+o*a)-l}function Wn(n,t,e,i,r){const s=Math.max(1e-6,e-t);if(n.length===0||!Number.isFinite(s))return{points:[],coverageRatio:0,centerCoverageRatio:0};const o=Math.max(1.5,Math.min(4,s/18)),a=Math.max(1,Math.ceil(s/o)),l=new Map;for(const f of n){const p=i(f);if(!Number.isFinite(p)||p<t||p>e)continue;const g=Math.max(0,Math.min(a-1,Math.floor((p-t)/o))),x=f.weight/(1+Math.abs(r(f))),y=l.get(g);(!y||x>y.score)&&l.set(g,{point:f,score:x})}const u=Array.from(l.values()).sort((f,p)=>i(f.point)-i(p.point)).map(f=>f.point),d=Math.max(0,Math.floor(a*.3)),c=Math.max(d+1,Math.ceil(a*.7));let m=0;for(let f=d;f<c;f++)l.has(f)&&m++;const h=Math.max(1,c-d);return{points:u,coverageRatio:u.length/a,centerCoverageRatio:m/h}}function jn(n,t){const e=n.dirX*t.dirY-n.dirY*t.dirX;if(!Number.isFinite(e)||Math.abs(e)<=1e-6)return null;const i=t.pointX-n.pointX,r=t.pointY-n.pointY,s=(i*t.dirY-r*t.dirX)/e;return{x:n.pointX+n.dirX*s,y:n.pointY+n.dirY*s}}function qo(n){if(n.length<3)return 0;let t=0;for(let e=0;e<n.length;e++){const i=n[e],r=n[(e+1)%n.length];t+=i.x*r.y-r.x*i.y}return t*.5}function Ko(n,t,e,i,r,s,o,a,l){const u=jo(i,r.label,s,o,r,t),d=u.map(N=>({x:N.x,y:N.y}));if(u.length<l.minEdgePoints)return{candidate:null,failureStage:"min_edge_points",pointsCount:u.length,strongEdgePoints:d};const c=Ho(u);if(!c)return{candidate:null,failureStage:"dominant_axes",pointsCount:u.length,strongEdgePoints:d};const m=u.map(N=>{const lt=N.x-c.centerX,ht=N.y-c.centerY;return{...N,u:lt*c.dirX+ht*c.dirY,v:lt*c.orthoX+ht*c.orthoY}}),h={x:c.centerX,y:c.centerY},f=Te(m.map(N=>({value:N.u,weight:N.weight})),l.extentQuantileLow),p=Te(m.map(N=>({value:N.u,weight:N.weight})),l.extentQuantileHigh),g=Te(m.map(N=>({value:N.v,weight:N.weight})),l.extentQuantileLow),x=Te(m.map(N=>({value:N.v,weight:N.weight})),l.extentQuantileHigh),y=Math.max(1e-6,Math.max(Math.abs(f),Math.abs(p))),b=Math.max(1e-6,Math.max(Math.abs(g),Math.abs(x))),_=72,F=360/_,v=Array.from({length:_},()=>[]),w=N=>{let lt=N%360;return lt<0&&(lt+=360),lt},P=(N,lt)=>{const ht=Math.abs(w(N)-w(lt));return Math.min(ht,360-ht)};m.forEach(N=>{const lt=N.u/y,ht=N.v/b,At=w(Math.atan2(ht,lt)*180/Math.PI),Yt=Math.hypot(lt,ht),Nt=Math.max(0,Math.min(_-1,Math.floor(At/F)));v[Nt].push({point:N,angleDeg:At,normRadius:Yt})});const T=v.map(N=>N.length>0?Br(N.map(lt=>lt.normRadius)):-1/0),C=(N,lt)=>{let ht=-1,At=-1/0;for(let yt=0;yt<v.length;yt++){if(v[yt].length===0)continue;const me=(yt+.5)*F;if(P(me,N)>45||lt.some(Bs=>P(me,Bs)<45))continue;const Ae=T[yt];Ae>At&&(At=Ae,ht=yt)}let Yt=ht>=0?(ht+.5)*F:N,Nt=ht>=0?v[ht]:m.map(yt=>{const pe=yt.u/y,me=yt.v/b;return{point:yt,angleDeg:w(Math.atan2(me,pe)*180/Math.PI),normRadius:Math.hypot(pe,me)}}).filter(yt=>P(yt.angleDeg,N)<=45&&!lt.some(pe=>P(yt.angleDeg,pe)<45));if(Nt.length===0&&(Nt=m.map(yt=>{const pe=yt.u/y,me=yt.v/b;return{point:yt,angleDeg:w(Math.atan2(me,pe)*180/Math.PI),normRadius:Math.hypot(pe,me)}}).filter(yt=>P(yt.angleDeg,N)<=45),Yt=N),Nt.length===0)return{x:m[0].x,y:m[0].y,u:m[0].u,v:m[0].v,angleDeg:N};const On=ht>=0?T[ht]:Br(Nt.map(yt=>yt.normRadius));let we=0,Mn=0,hr=0,dr=0,fr=0;for(const yt of Nt){const pe=P(yt.angleDeg,N)/45,me=Math.abs(yt.normRadius-On),Ae=Math.max(1e-6,yt.point.weight)/(1+pe*2+me*6);we+=Ae,Mn+=yt.point.x*Ae,hr+=yt.point.y*Ae,dr+=yt.point.u*Ae,fr+=yt.point.v*Ae}return we>0?{x:Mn/we,y:hr/we,u:dr/we,v:fr/we,angleDeg:Yt}:{x:Nt[0].point.x,y:Nt[0].point.y,u:Nt[0].point.u,v:Nt[0].point.v,angleDeg:Nt[0].angleDeg}},k=C(225,[]),M=C(315,[k.angleDeg]),S=C(45,[k.angleDeg,M.angleDeg]),R=C(135,[k.angleDeg,M.angleDeg,S.angleDeg]),A=[{x:k.x,y:k.y},{x:M.x,y:M.y},{x:S.x,y:S.y},{x:R.x,y:R.y}],L=p-f,X=x-g,D=Math.min(L,X),E=Math.max(L,X);if(!Number.isFinite(D)||D<l.minSpanPx||E/Math.max(1,D)>l.maxAspectRatio)return{candidate:null,failureStage:"span_aspect",pointsCount:u.length,minSpan:D,maxSpan:E,axisCentroid:h,axisExtremePoints:A,strongEdgePoints:d};const U=Math.max(l.bandMinPx,Math.min(l.bandMaxPx,D*l.bandScale)),O=Math.max(1,Math.min(3,U*.55)),z=Math.max(a,0),J=void 0,V=void 0,Q=N=>N.map(lt=>({x:lt.x,y:lt.y,weight:lt.weight})),W=N=>N.map(lt=>({x:lt.x,y:lt.y})),tt=(N,lt,ht)=>N.filter(At=>{if(!Number.isFinite(At.weight)||At.weight<z)return!1;const Yt=Qo(n,t,e,At.x,At.y,lt,ht,O);return Number.isFinite(Yt)&&Yt>=l.minPointContrast}),at=f,B=p,Z=g,q=x,I=l.minCoverageRatio,H=l.minCenterCoverageRatio,G=[],Y=[],et=[],nt=[],it=[],ct=(N,lt,ht,At,Yt,Nt)=>(ht-N)*(Nt-lt)-(At-lt)*(Yt-N),dt=N=>N>1e-6?1:N<-1e-6?-1:0,rt=[{u:(k.u+M.u)*.5,v:(k.v+M.v)*.5},{u:(M.u+S.u)*.5,v:(M.v+S.v)*.5},{u:(S.u+R.u)*.5,v:(S.v+R.v)*.5},{u:(R.u+k.u)*.5,v:(R.v+k.v)*.5}],Dt=(N,lt)=>{const ht=dt(ct(k.u,k.v,S.u,S.v,N,lt)),At=dt(ct(M.u,M.v,R.u,R.v,N,lt));return`${ht},${At}`},Gt=new Map;rt.forEach((N,lt)=>{Gt.set(Dt(N.u,N.v),lt)});for(const N of m){if(!Number.isFinite(N.u)||!Number.isFinite(N.v)){it.push(N);continue}let ht=Gt.get(Dt(N.u,N.v))??-1;if(ht<0){let At=1/0;for(let Yt=0;Yt<rt.length;Yt++){const Nt=rt[Yt],On=(N.u-Nt.u)/y,we=(N.v-Nt.v)/b,Mn=On*On+we*we;Mn<At&&(At=Mn,ht=Yt)}}ht===0?G.push(N):ht===1?Y.push(N):ht===2?et.push(N):ht===3?nt.push(N):it.push(N)}const j=[...G,...et],ut=[...Y,...nt],ot={dir:j.length,ortho:ut.length,unassigned:m.length-j.length-ut.length},pt=G.length>=l.minSidePoints?Te(G.map(N=>({value:N.v,weight:N.weight})),.5):g,mt=et.length>=l.minSidePoints?Te(et.map(N=>({value:N.v,weight:N.weight})),.5):x,ee=nt.length>=l.minSidePoints?Te(nt.map(N=>({value:N.u,weight:N.weight})),.5):f,Bt=Y.length>=l.minSidePoints?Te(Y.map(N=>({value:N.u,weight:N.weight})),.5):p,_t=[{x:(k.x+M.x)*.5,y:(k.y+M.y)*.5},{x:(M.x+S.x)*.5,y:(M.y+S.y)*.5},{x:(S.x+R.x)*.5,y:(S.y+R.y)*.5},{x:(R.x+k.x)*.5,y:(R.y+k.y)*.5}],ce=G.filter(N=>Math.abs(N.v-pt)<=U&&N.u>=at&&N.u<=B),$t=et.filter(N=>Math.abs(N.v-mt)<=U&&N.u>=at&&N.u<=B),Fe=nt.filter(N=>Math.abs(N.u-ee)<=U&&N.v>=Z&&N.v<=q),$=Y.filter(N=>Math.abs(N.u-Bt)<=U&&N.v>=Z&&N.v<=q),Ft=[ce.length,$.length,$t.length,Fe.length],Tt=[W(ce),W($),W($t),W(Fe)],xt=tt(ce,-c.orthoX,-c.orthoY),wt=tt($t,c.orthoX,c.orthoY),ne=tt(Fe,-c.dirX,-c.dirY),Jt=tt($,c.dirX,c.dirY),bn=[xt.length,Jt.length,wt.length,ne.length],Be=[W(xt),W(Jt),W(wt),W(ne)],Oe=Wn(xt,at,B,N=>N.u,N=>N.v-pt),Ve=Wn(Jt,Z,q,N=>N.v,N=>N.u-Bt),Ge=Wn(wt,at,B,N=>N.u,N=>N.v-mt),Xe=Wn(ne,Z,q,N=>N.v,N=>N.u-ee),An=(N,lt)=>N.slice().sort((ht,At)=>lt(ht)-lt(At)),Tn=An(xt,N=>N.u),In=An(Jt,N=>N.v),Rn=An(wt,N=>N.u),Nn=An(ne,N=>N.v),hi=[Tn.length,In.length,Rn.length,Nn.length],Ht=[Oe.coverageRatio,Ve.coverageRatio,Ge.coverageRatio,Xe.coverageRatio];Oe.centerCoverageRatio,Ve.centerCoverageRatio,Ge.centerCoverageRatio,Xe.centerCoverageRatio;const Xt=[W(Tn),W(In),W(Rn),W(Nn)],Zt={axisPointCounts:ot,sideBandPointCounts:Ft,sideContrastPointCounts:bn,gradientThreshold:a,pointAxisMinDot:J,pointAxisMargin:V,bandWidth:U,minPointContrast:l.minPointContrast,minCoverageRatio:I,minCenterCoverageRatio:H,axisCentroid:h,axisExtremePoints:A,axisSideCenters:_t,strongEdgePoints:d,axisDirPoints:W(j),axisOrthoPoints:W(ut),axisUnassignedPoints:W(it),sideBandPoints:Tt,sideContrastPoints:Be};if(Tn.length<l.minSidePoints||Rn.length<l.minSidePoints||Nn.length<l.minSidePoints||In.length<l.minSidePoints)return{candidate:null,failureStage:"min_side_points",pointsCount:u.length,minSpan:D,maxSpan:E,sidePointCounts:hi,sideCoverageRatios:Ht,...Zt,sideFitPoints:Xt};if(Oe.coverageRatio<I||Ve.coverageRatio<I||Ge.coverageRatio<I||Xe.coverageRatio<I||Oe.centerCoverageRatio<H||Ve.centerCoverageRatio<H||Ge.centerCoverageRatio<H||Xe.centerCoverageRatio<H)return{candidate:null,failureStage:"side_coverage",pointsCount:u.length,minSpan:D,maxSpan:E,sidePointCounts:hi,sideCoverageRatios:Ht,...Zt,sideFitPoints:Xt};const di=Tn,fi=Rn,pi=Nn,mi=In,ie=hi,ze=Pe(Q(di)),Ye=Pe(Q(fi)),We=Pe(Q(pi)),je=Pe(Q(mi));if(!ze||!Ye||!We||!je)return{candidate:null,failureStage:"fit_lines",pointsCount:u.length,minSpan:D,maxSpan:E,sidePointCounts:ie,sideCoverageRatios:Ht,...Zt,sideFitPoints:Xt};const Me=Co(Q(di),ze,Q(mi),je,Q(fi),Ye,Q(pi),We),Ln=l.minAxisDot,En=(N,lt,ht)=>Math.abs(N.dirX*lt+N.dirY*ht),zt=[En(ze,c.dirX,c.dirY),En(je,c.orthoX,c.orthoY),En(Ye,c.dirX,c.dirY),En(We,c.orthoX,c.orthoY)];if(zt[0]<Ln||zt[1]<Ln||zt[2]<Ln||zt[3]<Ln)return{candidate:null,failureStage:"axis_alignment",pointsCount:u.length,minSpan:D,maxSpan:E,sidePointCounts:ie,sideCoverageRatios:Ht,axisDots:zt,...Zt,sideFitPoints:Xt,sideFitLines:Me};const re=Math.max(l.residualLimitFloor,U*l.residualLimitScale),gt=[Yn(Q(di),ze),Yn(Q(fi),Ye),Yn(Q(pi),We),Yn(Q(mi),je)],gi=[gt[0],gt[3],gt[1],gt[2]],yi=Math.max(...gt),ue=jn(ze,We),he=jn(ze,je),de=jn(Ye,je),fe=jn(Ye,We);if(!ue||!he||!de||!fe)return{candidate:null,failureStage:"corners",pointsCount:u.length,minSpan:D,maxSpan:E,sidePointCounts:ie,sideCoverageRatios:Ht,axisDots:zt,sideResiduals:[gt[0],gt[3],gt[1],gt[2]],residualLimit:re,...Zt,sideFitPoints:Xt,sideFitLines:Me};const _n=[ue,he,de,fe],se=Math.abs(qo(_n));if(!Number.isFinite(se)||se<l.minQuadArea)return{candidate:null,failureStage:"quad_area",pointsCount:u.length,minSpan:D,maxSpan:E,sidePointCounts:ie,sideCoverageRatios:Ht,axisDots:zt,sideResiduals:[gt[0],gt[3],gt[1],gt[2]],residualLimit:re,...Zt,sideFitPoints:Xt,sideFitLines:Me,quadArea:se};const xi=Math.hypot(he.x-ue.x,he.y-ue.y),bi=Math.hypot(de.x-he.x,de.y-he.y),_i=Math.hypot(de.x-fe.x,de.y-fe.y),Mi=Math.hypot(fe.x-ue.x,fe.y-ue.y),He=[xi,bi,_i,Mi],wi=Math.min(xi,bi,_i,Mi),Ls=Math.max(xi,bi,_i,Mi);if(!Number.isFinite(wi)||wi<l.minSideLength||Ls/Math.max(1,wi)>l.maxAspectRatio)return{candidate:null,failureStage:"side_length",pointsCount:u.length,minSpan:D,maxSpan:E,sidePointCounts:ie,sideCoverageRatios:Ht,axisDots:zt,sideResiduals:[gt[0],gt[3],gt[1],gt[2]],residualLimit:re,...Zt,sideFitPoints:Xt,sideFitLines:Me,quadArea:se,sideLengths:He};const Qe=ks(he.x-ue.x+(de.x-fe.x),he.y-ue.y+(de.y-fe.y));if(!Qe)return{candidate:null,failureStage:"corners",pointsCount:u.length,minSpan:D,maxSpan:E,sidePointCounts:ie,sideCoverageRatios:Ht,axisDots:zt,sideResiduals:[gt[0],gt[3],gt[1],gt[2]],residualLimit:re,...Zt,sideFitPoints:Xt,sideFitLines:Me,quadArea:se,sideLengths:He};const lr={x:-Qe.y,y:Qe.x},Si=(ue.x+he.x+de.x+fe.x)*.25,vi=(ue.y+he.y+de.y+fe.y)*.25,Un=_n.map(N=>{const lt=N.x-Si,ht=N.y-vi;return{u:lt*Qe.x+ht*Qe.y,v:lt*lr.x+ht*lr.y}}),Dn=(Math.max(...Un.map(N=>N.u))-Math.min(...Un.map(N=>N.u)))*.5,Bn=(Math.max(...Un.map(N=>N.v))-Math.min(...Un.map(N=>N.v)))*.5;if(!Number.isFinite(Dn)||!Number.isFinite(Bn)||Math.min(Dn,Bn)<6)return{candidate:null,failureStage:"box_size",pointsCount:u.length,minSpan:D,maxSpan:E,sidePointCounts:ie,sideCoverageRatios:Ht,axisDots:zt,sideResiduals:[gt[0],gt[3],gt[1],gt[2]],residualLimit:re,sideFitPoints:Xt,quadArea:se,sideLengths:He};const K=il(n,t,e,_n,ze,je,Ye,We,Si,vi,Dn,Bn,l.innerPurityStdScale,l.outerMeanSpreadLimit);if(!Number.isFinite(yi)||yi>re)return{candidate:null,failureStage:"residual",pointsCount:u.length,minSpan:D,maxSpan:E,sidePointCounts:ie,sideCoverageRatios:Ht,axisDots:zt,sideResiduals:gi,residualLimit:re,...Zt,sideFitPoints:Xt,sideFitLines:Me,quadArea:se,sideLengths:He,outerContrast:K.contrast,outerUniformityOk:K.ok,outerMeanSpread:K.meanSpread,outerMeanSpreadLimit:K.meanSpreadLimit,outerAvgStd:K.avgStd,outerAvgStdLimit:K.avgStdLimit,outerSideMeans:K.outerSideMeans,outerSideStds:K.outerSideStds,outerSideStdLimit:K.outerSideStdLimit,outerSideQuads:K.outerSideQuads,innerSideUniformityOk:K.innerSideOk,innerSideStds:K.innerSideStds,innerSideStdLimit:K.innerSideStdLimit,innerSideQuads:K.innerSideQuads};const cr=l.filterBlockPurity&&(!K.ok||!K.innerSideOk);if(cr||K.contrast<l.minOuterContrast)return{candidate:null,failureStage:cr?K.ok?"inner_roi_uniformity":"outer_uniformity":"outer_contrast",pointsCount:u.length,minSpan:D,maxSpan:E,sidePointCounts:ie,sideCoverageRatios:Ht,axisDots:zt,sideResiduals:gi,residualLimit:re,...Zt,sideFitPoints:Xt,sideFitLines:Me,quadArea:se,sideLengths:He,outerContrast:K.contrast,outerUniformityOk:K.ok,outerMeanSpread:K.meanSpread,outerMeanSpreadLimit:K.meanSpreadLimit,outerAvgStd:K.avgStd,outerAvgStdLimit:K.avgStdLimit,outerSideMeans:K.outerSideMeans,outerSideStds:K.outerSideStds,outerSideStdLimit:K.outerSideStdLimit,outerSideQuads:K.outerSideQuads,innerSideUniformityOk:K.innerSideOk,innerSideStds:K.innerSideStds,innerSideStdLimit:K.innerSideStdLimit,innerSideQuads:K.innerSideQuads};const ur=Mt(Vt(_n,1),t,e);if(!ur)return{candidate:null,failureStage:"bbox",pointsCount:u.length,minSpan:D,maxSpan:E,sidePointCounts:ie,sideCoverageRatios:Ht,axisDots:zt,sideResiduals:[gt[0],gt[3],gt[1],gt[2]],residualLimit:re,...Zt,sideFitPoints:Xt,sideFitLines:Me,quadArea:se,sideLengths:He,outerContrast:K.contrast,outerUniformityOk:K.ok,outerMeanSpread:K.meanSpread,outerMeanSpreadLimit:K.meanSpreadLimit,outerAvgStd:K.avgStd,outerAvgStdLimit:K.avgStdLimit,outerSideMeans:K.outerSideMeans,outerSideStds:K.outerSideStds,outerSideStdLimit:K.outerSideStdLimit,outerSideQuads:K.outerSideQuads,innerSideUniformityOk:K.innerSideOk,innerSideStds:K.innerSideStds,innerSideStdLimit:K.innerSideStdLimit,innerSideQuads:K.innerSideQuads};const Es=1/(1+yi/Math.max(1,re)),Us=l.filterBlockPurity?K.score:1,Ds=K.contrast*Us*Es*Math.sqrt(se);return{candidate:{centerX:Si,centerY:vi,dirX:Qe.x,dirY:Qe.y,halfWidth:Dn,halfHeight:Bn,score:Ds,bbox:ur,corners:_n,sideFitPoints:Xt,outerSideMeans:K.outerSideMeans,outerSideQuads:K.outerSideQuads},failureStage:null,pointsCount:u.length,minSpan:D,maxSpan:E,sidePointCounts:ie,sideCoverageRatios:Ht,axisDots:zt,sideResiduals:gi,residualLimit:re,...Zt,sideFitPoints:Xt,sideFitLines:Me,quadArea:se,sideLengths:He,outerContrast:K.contrast,outerUniformityOk:K.ok,outerMeanSpread:K.meanSpread,outerMeanSpreadLimit:K.meanSpreadLimit,outerAvgStd:K.avgStd,outerAvgStdLimit:K.avgStdLimit,outerSideMeans:K.outerSideMeans,outerSideStds:K.outerSideStds,outerSideStdLimit:K.outerSideStdLimit,outerSideQuads:K.outerSideQuads,innerSideUniformityOk:K.innerSideOk,innerSideStds:K.innerSideStds,innerSideStdLimit:K.innerSideStdLimit,innerSideQuads:K.innerSideQuads}}function $o(n,t,e,i,r,s,o,a,l){return Ko(n,t,e,i,r,s,o,a,l).candidate}function Jo(n,t,e,i,r,s){const o=nr(r),a=Math.max(i*8,i+128),l=tl(n,t,e,o.downsampleMaxSide);s==null||s("Detecting candidates: downsampling...",.02),s==null||s("Detecting candidates: edge stage...",.06);const u=zo(l.gray,l.width,l.height),d=u.gray,c=u.gradient;s==null||s(`Detecting candidates: gradient (${u.backend==="webgl"?"WebGL1":"CPU"})...`,.1);const m=Vo(c.magnitude,o.gradientPercentiles),h=l.width*l.height,f=Math.max(o.minComponentAreaPx,Math.round(h*o.minComponentAreaRatio)),p=Math.max(f+1,Math.round(h*o.maxComponentAreaRatio)),g=[],x=Math.max(1,m.reduce((k,M,S)=>k+(S<=1,2),0));let y=0;for(let k=0;k<m.length;k++){const M=m[k],S=Go(c.magnitude,M),R=k<=1?[3,2]:[2,1];for(const A of R){const L=y/x;s==null||s(`Detecting candidates: threshold ${k+1}/${m.length}, dilate ${A}`,.12+.78*L);const X=Yo(S,l.width,l.height,A),{labels:D,components:E}=Wo(X,l.width,l.height);for(const U of E){if(U.touchesBorder||U.area<f||U.area>p)continue;const O=$o(d,l.width,l.height,D,U,S,c,M,o);if(!O)continue;const z=1/l.scale,J=O.corners.map(V=>({x:V.x*z,y:V.y*z}));g.push({centerX:O.centerX*z,centerY:O.centerY*z,dirX:O.dirX,dirY:O.dirY,halfWidth:O.halfWidth*z,halfHeight:O.halfHeight*z,score:O.score,bbox:{x:O.bbox.x*z,y:O.bbox.y*z,w:O.bbox.w*z,h:O.bbox.h*z},corners:J,sideFitPoints:O.sideFitPoints?[O.sideFitPoints[0].map(V=>({x:V.x*z,y:V.y*z})),O.sideFitPoints[1].map(V=>({x:V.x*z,y:V.y*z})),O.sideFitPoints[2].map(V=>({x:V.x*z,y:V.y*z})),O.sideFitPoints[3].map(V=>({x:V.x*z,y:V.y*z}))]:void 0,outerSideMeans:O.outerSideMeans,outerSideQuads:O.outerSideQuads?[O.outerSideQuads[0].map(V=>({x:V.x*z,y:V.y*z})),O.outerSideQuads[1].map(V=>({x:V.x*z,y:V.y*z})),O.outerSideQuads[2].map(V=>({x:V.x*z,y:V.y*z})),O.outerSideQuads[3].map(V=>({x:V.x*z,y:V.y*z}))]:void 0})}g.length>a&&(g.sort((U,O)=>O.score-U.score),g.length=a),y++}}console.log(`[SFR Auto Detect] Candidate pool before dedupe: ${g.length}`),s==null||s(`Detecting candidates: deduplicating (0/${Math.max(1,Math.min(g.length,Math.max(i*4,i+32)))})...`,.94),g.sort((k,M)=>M.score-k.score);const b=Math.max(i*4,i+32),_=g.length>b?g.slice(0,b):g,F=[];if(_.length<=256){console.log(`[SFR Auto Detect] Using simple dedupe for ${_.length} candidates`);for(let k=0;k<_.length;k++){const M=_[k];console.log(`[SFR Auto Detect] Simple dedupe candidate ${k+1}/${_.length}`,M.bbox);const S=_.length<=0?1:k/_.length;if(s==null||s(`Detecting candidates: deduplicating (${k}/${_.length})...`,.94+.05*Math.min(1,S)),!F.some(A=>{const L=Math.hypot(M.centerX-A.centerX,M.centerY-A.centerY),X=Math.max(Math.hypot(M.bbox.w,M.bbox.h),Math.hypot(A.bbox.w,A.bbox.h));return Vr(M.bbox,A.bbox)>.28||L<X*.18})&&(F.push(M),F.length>=i))break}return s==null||s("Detecting candidates: deduplicating...",1),F}const v=Math.max(32,Math.round(Math.sqrt(Math.max(1,t*e)/4096))),w=new Map,P=new Set,T=k=>Math.floor(k/v),C=(k,M)=>{if(!Number.isFinite(k.bbox.x)||!Number.isFinite(k.bbox.y)||!Number.isFinite(k.bbox.w)||!Number.isFinite(k.bbox.h)||k.bbox.w<=0||k.bbox.h<=0||k.bbox.w>t*4||k.bbox.h>e*4)return;const S=T(k.bbox.x),R=T(k.bbox.x+k.bbox.w),A=T(k.bbox.y),L=T(k.bbox.y+k.bbox.h);for(let X=A;X<=L;X++)for(let D=S;D<=R;D++){const E=`${D},${X}`,U=w.get(E);U?U.push(M):w.set(E,[M])}};for(let k=0;k<_.length;k++){const M=_[k];if(k===0||k%200===0){const D=_.length<=0?1:k/_.length;s==null||s(`Detecting candidates: deduplicating (${k}/${_.length})...`,.94+.05*Math.min(1,D))}P.clear();const S=T(M.bbox.x),R=T(M.bbox.x+M.bbox.w),A=T(M.bbox.y),L=T(M.bbox.y+M.bbox.h);let X=!1;for(let D=A-1;D<=L+1&&!X;D++)for(let E=S-1;E<=R+1&&!X;E++){const U=w.get(`${E},${D}`);if(U)for(const O of U){if(P.has(O))continue;P.add(O);const z=F[O];if(!z)continue;const J=Math.hypot(M.centerX-z.centerX,M.centerY-z.centerY),V=Math.max(Math.hypot(M.bbox.w,M.bbox.h),Math.hypot(z.bbox.w,z.bbox.h));if(Vr(M.bbox,z.bbox)>.28||J<V*.18){X=!0;break}}}if(!X){const D=F.length;if(F.push(M),C(M,D),F.length>=i)break}}return s==null||s("Detecting candidates: deduplicating...",1),F}function Zo(n,t,e){const i=new Uint8Array(n.length);for(let r=0;r<e;r++)for(let s=0;s<t;s++){let o=0,a=0;for(let l=-1;l<=1;l++){const u=r+l;if(!(u<0||u>=e))for(let d=-1;d<=1;d++){const c=s+d;c<0||c>=t||(o+=n[u*t+c],a++)}}i[r*t+s]=Math.round(o/Math.max(1,a))}return i}function tl(n,t,e,i){const r=Math.max(t,e);if(r<=i)return{gray:n,width:t,height:e,scale:1};const s=i/r,o=Math.max(1,Math.round(t*s)),a=Math.max(1,Math.round(e*s)),l=new Uint8Array(o*a);for(let u=0;u<a;u++){const d=Math.min(e-1,Math.floor(u/s));for(let c=0;c<o;c++){const m=Math.min(t-1,Math.floor(c/s));l[u*o+c]=n[d*t+m]}}return{gray:l,width:o,height:a,scale:s}}function Vr(n,t){const e=Math.max(n.x,t.x),i=Math.max(n.y,t.y),r=Math.min(n.x+n.w,t.x+t.w),s=Math.min(n.y+n.h,t.y+t.h),o=Math.max(0,r-e),a=Math.max(0,s-i),l=o*a;if(l<=0)return 0;const u=n.w*n.h+t.w*t.h-l;return u>0?l/u:0}function Gr(n){const t=n.length;if(t===0)return{count:0,mean:0,std:1/0};let e=0;for(const s of n)e+=s;const i=e/t;let r=0;for(const s of n){const o=s-i;r+=o*o}return r/=t,{count:t,mean:i,std:Math.sqrt(Math.max(0,r))}}function el(n,t,e,i){return{p1:{x:n.x-t*i,y:n.y-e*i},p2:{x:n.x+t*i,y:n.y+e*i}}}function Xr(n,t,e,i,r){return[{x:n.p1.x+t*i,y:n.p1.y+e*i},{x:n.p2.x+t*i,y:n.p2.y+e*i},{x:n.p2.x+t*r,y:n.p2.y+e*r},{x:n.p1.x+t*r,y:n.p1.y+e*r}]}function nl(n,t,e){let i=0;for(let r=0;r<4;r++){const s=e[r],o=e[(r+1)%4],a=(o.x-s.x)*(t-s.y)-(o.y-s.y)*(n-s.x);if(Math.abs(a)<=1e-6)continue;const l=a>0?1:-1;if(i===0)i=l;else if(i!==l)return!1}return!0}function zr(n,t,e,i){const r=Mt(Vt(i,1),t,e);if(!r)return[];const s=[];for(let o=r.y;o<r.y+r.h;o++)for(let a=r.x;a<r.x+r.w;a++)nl(a,o,i)&&s.push(n[o*t+a]);return s}function il(n,t,e,i,r,s,o,a,l,u,d,c,m,h){const f=d*2,p=c*2,g=Math.hypot(i[1].x-i[0].x,i[1].y-i[0].y),x=Math.hypot(i[2].x-i[1].x,i[2].y-i[1].y),y=Math.hypot(i[2].x-i[3].x,i[2].y-i[3].y),b=Math.hypot(i[3].x-i[0].x,i[3].y-i[0].y),F=Math.max(...[g,x,y,b]),v=Math.max(2,Math.min(f,p)),w=Math.max(4,F*.25),P=Math.max(2,Math.min(12,v*.22)),T=Math.max(1,Math.min(P,Math.max(1,v*.5-1))),C=1,k=Math.max(8,Math.round(Math.min(w,P*3))),M=[[i[0],i[1],i[1],i[0]],[i[1],i[2],i[2],i[1]],[i[2],i[3],i[3],i[2]],[i[3],i[0],i[0],i[3]]],S=[[i[0],i[1],i[1],i[0]],[i[1],i[2],i[2],i[1]],[i[2],i[3],i[3],i[2]],[i[3],i[0],i[0],i[3]]],R=[],A=[],L=[{corners:[i[0],i[1]],seedLine:r,sideLength:g},{corners:[i[1],i[2]],seedLine:s,sideLength:x},{corners:[i[2],i[3]],seedLine:o,sideLength:y},{corners:[i[3],i[0]],seedLine:a,sideLength:b}];for(let G=0;G<L.length;G++){const Y=L[G],et=Math.max(1,Y.sideLength*.5-1),nt=Math.max(1,Math.min(et,w*.5)),it={x:(Y.corners[0].x+Y.corners[1].x)*.5,y:(Y.corners[0].y+Y.corners[1].y)*.5},ct=el(it,Y.seedLine.dirX,Y.seedLine.dirY,nt),dt=ye(n,t,e,ct.p1,ct.p2,nt,Math.max(4,C+Math.max(P,T)+2)),rt=(dt==null?void 0:dt.line)||ct,Dt=rt.p2.x-rt.p1.x,Gt=rt.p2.y-rt.p1.y,j=Math.hypot(Dt,Gt);if(!Number.isFinite(j)||j<=1e-6)return{ok:!1,score:0,meanSpread:1/0,meanSpreadLimit:1/0,avgStd:1/0,avgStdLimit:1/0,contrast:0,outerMean:0,outerSideMeans:[0,0,0,0],outerSideStds:[1/0,1/0,1/0,1/0],outerSideStdLimit:1/0,outerSideQuads:M,innerSideOk:!1,innerSideStds:[1/0,1/0,1/0,1/0],innerSideStdLimit:1/0,innerSideQuads:S};let ut=-Gt/j,ot=Dt/j;const pt={x:(rt.p1.x+rt.p2.x)*.5,y:(rt.p1.y+rt.p2.y)*.5};(pt.x-l)*ut+(pt.y-u)*ot<0&&(ut=-ut,ot=-ot);const mt=Xr(rt,ut,ot,C,C+P),ee=Xr(rt,ut,ot,-C,-(C+T));M[G]=mt,S[G]=ee,R.push(zr(n,t,e,mt)),A.push(zr(n,t,e,ee))}const X=R.map(Gr);if(X.some(G=>G.count<k||!Number.isFinite(G.std)))return{ok:!1,score:0,meanSpread:1/0,meanSpreadLimit:1/0,avgStd:1/0,avgStdLimit:1/0,contrast:0,outerMean:0,outerSideMeans:[0,0,0,0],outerSideStds:[1/0,1/0,1/0,1/0],outerSideStdLimit:1/0,outerSideQuads:M,innerSideOk:!1,innerSideStds:[1/0,1/0,1/0,1/0],innerSideStdLimit:1/0,innerSideQuads:S};const D=A.map(Gr);if(D.some(G=>G.count<k||!Number.isFinite(G.std)||!Number.isFinite(G.mean)))return{ok:!1,score:0,meanSpread:1/0,meanSpreadLimit:1/0,avgStd:1/0,avgStdLimit:1/0,contrast:0,outerMean:0,outerSideMeans:[0,0,0,0],outerSideStds:[1/0,1/0,1/0,1/0],outerSideStdLimit:1/0,outerSideQuads:M,innerSideOk:!1,innerSideStds:[1/0,1/0,1/0,1/0],innerSideStdLimit:1/0,innerSideQuads:S};const E=X.map(G=>G.mean),U=E.reduce((G,Y)=>G+Y,0)/E.length,O=D.reduce((G,Y)=>G+Y.mean,0)/D.length,z=Math.abs(O-U),J=Math.max(...E)-Math.min(...E),V=X.reduce((G,Y)=>G+Y.std,0)/X.length,Q=Math.max(0,h),W=Math.max(6,Math.min(20,z*.45)),tt=E,at=X.map(G=>G.std),B=Math.max(W,Math.min(30,W*m)),Z=D.map(G=>G.std),q=Z.every(G=>G<=B),I=J<=Q&&V<=W,H=1/(1+J/Math.max(1,Q)+V/Math.max(1,W));return{ok:I,score:H,meanSpread:J,meanSpreadLimit:Q,avgStd:V,avgStdLimit:W,contrast:z,outerMean:U,outerSideMeans:tt,outerSideStds:at,outerSideStdLimit:W,outerSideQuads:M,innerSideOk:q,innerSideStds:Z,innerSideStdLimit:B,innerSideQuads:S}}function ye(n,t,e,i,r,s,o){const a=r.x-i.x,l=r.y-i.y,u=Math.hypot(a,l);if(!Number.isFinite(u)||u<=1e-6)return null;const d=a/u,c=l/u,m=-c,h=d,f=(i.x+r.x)*.5,p=(i.y+r.y)*.5,g=be({p1:i,p2:r},o+2);if(!Mt(Vt(g||[i,r],2),t,e))return null;const y=Math.max(8,Math.round(s*2)+1),b=Math.max(8,Math.round(o*2)+1),_=y>1?s*2/(y-1):0,F=b>1?o*2/(b-1):0,v=Array.from({length:y},()=>new Array(b).fill(0));for(let M=0;M<y;M++){const S=-s+_*M;for(let R=0;R<b;R++){const A=-o+F*R,L=f+S*d+A*m,X=p+S*c+A*h;v[M][R]=rl(n,t,e,L,X)}}const w=vs(v,-o,-s,F,_,!0);if(w.length<8)return null;const P=w.map(M=>{const S=M.x,R=M.y;return{x:f+R*d+S*m,y:p+R*c+S*h,weight:M.weight}}),T=Pe(P);if(!T)return null;let C=T.dirX,k=T.dirY;return C*d+k*c<0&&(C=-C,k=-k),{line:{p1:{x:T.pointX-C*s,y:T.pointY-k*s},p2:{x:T.pointX+C*s,y:T.pointY+k*s}},fitPoints:P.map(M=>({x:M.x,y:M.y}))}}function rl(n,t,e,i,r){if(t<=0||e<=0||n.length!==t*e)return 0;const s=Math.max(0,Math.min(t-1,i)),o=Math.max(0,Math.min(e-1,r)),a=Math.floor(s),l=Math.floor(o),u=Math.min(t-1,a+1),d=Math.min(e-1,l+1),c=s-a,m=o-l,h=n[l*t+a],f=n[l*t+u],p=n[d*t+a],g=n[d*t+u],x=h*(1-c)+f*c,y=p*(1-c)+g*c;return x*(1-m)+y*m}function sl(n,t,e,i,r){if(i<=0||r<=0||i>=t-1||r>=e-1)return{gx:0,gy:0};const s=r*t+i;return{gx:(n[s+1]-n[s-1])*.5,gy:(n[s+t]-n[s-t])*.5}}function al(n){if(n.length<20)return null;const t=n.map(w=>Math.max(0,w.weight));let e=0;for(const w of t)e=Math.max(e,w);if(!(e>0))return null;for(let w=0;w<t.length;w++)t[w]/=e;const i=w=>{let P=0,T=0,C=0;for(let W=0;W<n.length;W++){const tt=w[W];tt>0&&(P+=tt,T+=n[W].x*tt,C+=n[W].y*tt)}if(!(P>0))return null;T/=P,C/=P;let k=0,M=0,S=0;for(let W=0;W<n.length;W++){const tt=w[W];if(!(tt>0))continue;const at=n[W].x-T,B=n[W].y-C;k+=tt*at*at,M+=tt*at*B,S+=tt*B*B}k/=P,M/=P,S/=P;const R=k+S,A=k*S-M*M,L=-R,X=A,D=Math.max(0,L*L-4*X),E=-.5*(L+(L>=0?1:-1)*Math.sqrt(D)),U=Math.abs(E)>1e-12?E:0,O=Math.abs(E)>1e-12?X/E:R,z=Math.max(U,O);let J=0,V=1;Math.abs(M)>1e-10?(J=z-S,V=M):k>S&&(J=1,V=0);const Q=Math.atan2(-J,V);return{centroid:{x:T,y:C},angle:Q,totalWeight:P}},r=i(t);if(!r)return null;const s=Math.cos(r.angle),o=Math.sin(r.angle),a=new Array(2*16*8).fill(0);for(let w=0;w<n.length;w++){const P=n[w].x-r.centroid.x,T=n[w].y-r.centroid.y,C=P*s+T*o,k=Math.round(C*8+16*8);if(k>=3&&k<a.length-3)for(let M=-3;M<=3;M++)a[k+M]+=t[w]}let l=16*8;for(let w=-5*8+16*8;w<=5*8+16*8;w++)a[w]>a[l]&&(l=w);let u=l-1;for(;u>1&&a[u]>.05*a[l];)u--;let d=l+1;for(;d<a.length-1&&a[d]>.05*a[l];)d++;let c=Math.max(1,u-8);for(;c>1&&a[c]<=a[u];)c--;let m=Math.min(a.length-1,d+8);for(;m<a.length-1&&a[m]<=a[d];)m++;const h=a.slice();for(let w=1;w<h.length;w++)h[w]+=h[w-1];const f=h[h.length-1];if(!(f>0))return null;let p=0;for(let w=1;w<h.length;w++)Math.abs(h[w]-.1*f)<Math.abs(h[p]-.1*f)&&(p=w);let g=h.length-1;for(let w=h.length-2;w>0;w--)Math.abs(h[w]-.9*f)<Math.abs(h[g]-.9*f)&&(g=w);let x=p/8-16,y=g/8-16;const b=y-x;x-=b*.7,y+=b*.7,x=Math.max((c+u)/16-16,x),y=Math.min((m+d)/16-16,y);const _=t.slice();for(let w=0;w<n.length;w++){const P=n[w].x-r.centroid.x,T=n[w].y-r.centroid.y,C=P*s+T*o;_[w]=C>=x&&C<=y?t[w]**4*(1/(10+Math.abs(C))):0}const F=i(_);if(!F)return null;const v=[];for(let w=0;w<n.length;w++)_[w]>0&&v.push({x:n[w].x,y:n[w].y,weight:_[w]});return v.length<8?null:{centroid:F.centroid,angle:F.angle,keptSamples:v}}function ol(n,t,e,i,r,s=Lt){var P;const o=r.x-i.x,a=r.y-i.y,l=Math.hypot(o,a);if(!Number.isFinite(l)||l<=12)return null;const u=o/l,d=a/l,c=5,m=4*s+.5,h=(T,C,k,M,S)=>{const R={x:T.x-C*l*.5,y:T.y-k*l*.5},A={p1:R,p2:{x:R.x+C*l,y:R.y+k*l}},L=be(A,m+2),X=Mt(Vt(L??[A.p1,A.p2],3),t,e),D=[],E=new Map;if(!X)return{reduced:null,scanlines:E};for(let U=X.y;U<X.y+X.h;U++)for(let O=X.x;O<X.x+X.w;O++){const z=O,J=U,V=z-R.x,Q=J-R.y,W=V*C+Q*k;if(!(W>c&&W<l-c))continue;const tt=z-T.x,at=J-T.y,B=tt*M+at*S;if(Math.abs(B)<12){const{gx:Z,gy:q}=sl(n,t,e,O,U),I=Z*Z+q*q;I>0&&D.push({x:z,y:J,weight:I})}if(Math.abs(B)<m){const Z=E.get(U);Z?(O<Z.start&&(Z.start=O),O>Z.end&&(Z.end=O)):E.set(U,{start:O,end:O})}}return{reduced:al(D),scanlines:E}};let f={x:(i.x+r.x)*.5,y:(i.y+r.y)*.5},p=u,g=d,x=-g,y=p,b=h(f,p,g,x,y);if(!b.reduced)return null;f=b.reduced.centroid,x=Math.cos(b.reduced.angle),y=Math.sin(b.reduced.angle),p=-y,g=x,p*u+g*d<0&&(p=-p,g=-g,x=-x,y=-y);let _=h(f,p,g,x,y);if(!_.reduced)return null;const F=Math.hypot(_.reduced.centroid.x-f.x,_.reduced.centroid.y-f.y);f=_.reduced.centroid,x=Math.cos(_.reduced.angle),y=Math.sin(_.reduced.angle),p=-y,g=x,p*u+g*d<0&&(p=-p,g=-g,x=-x,y=-y);let v=_;if(F>1){const T=h(f,p,g,x,y);T.reduced&&(v=T,f=T.reduced.centroid,x=Math.cos(T.reduced.angle),y=Math.sin(T.reduced.angle),p=-y,g=x,p*u+g*d<0&&(p=-p,g=-g))}const w=(((P=v.reduced)==null?void 0:P.keptSamples)??[]).map(T=>({x:T.x,y:T.y}));return w.length<8?null:{line:{p1:{x:f.x-p*l*.5,y:f.y-g*l*.5},p2:{x:f.x+p*l*.5,y:f.y+g*l*.5}},fitPoints:w,correctedScanlines:v.scanlines}}function ll(n,t,e){var b;const i=n.length,r=((b=n[0])==null?void 0:b.length)??0;if(r===0||i===0)return 0;const s=Math.max(0,Math.min(r-1,t)),o=Math.max(0,Math.min(i-1,e)),a=Math.floor(s),l=Math.floor(o),u=Math.min(r-1,a+1),d=Math.min(i-1,l+1),c=s-a,m=o-l,h=n[l][a],f=n[l][u],p=n[d][a],g=n[d][u],x=h*(1-c)+f*c,y=p*(1-c)+g*c;return x*(1-m)+y*m}function cl(n,t,e,i,r,s,o){var D;const a=i.p2.x-i.p1.x,l=i.p2.y-i.p1.y,u=Math.hypot(a,l);if(!Number.isFinite(u)||u<=1e-6)return null;const d=a/u,c=l/u,m=-c,h=d,f=(i.p1.x+i.p2.x)*.5,p=(i.p1.y+i.p2.y)*.5,g=be(i,s+2),x=Mt(Vt(g||[i.p1,i.p2],3),t,e);if(!x)return null;const y=bo(n,t,e,x,0,0,(o==null?void 0:o.bayerPattern)||"RGGB",o==null?void 0:o.greenPhase,o==null?void 0:o.blackLevel),b=y.length;if((((D=y[0])==null?void 0:D.length)??0)<6||b<6)return null;const F=Math.max(8,Math.round(r*2)+1),v=Math.max(8,Math.round(s*2)+1),w=F>1?r*2/(F-1):0,P=v>1?s*2/(v-1):0,T=Array.from({length:F},()=>new Array(v).fill(0));for(let E=0;E<F;E++){const U=-r+w*E;for(let O=0;O<v;O++){const z=-s+P*O,J=f+U*d+z*m,V=p+U*c+z*h;T[E][O]=ll(y,J-x.x,V-x.y)}}const{gx:C,gy:k}=_o(T),M=C>=k,S=vs(T,-s,-r,P,w,M);if(S.length<8)return null;const R=S.map(E=>{const U=E.x,O=E.y;return{x:f+O*d+U*m,y:p+O*c+U*h,weight:E.weight}}),A=Pe(R);if(!A)return null;let L=A.dirX,X=A.dirY;return L*d+X*c<0&&(L=-L,X=-X),{line:{p1:{x:A.pointX-L*r,y:A.pointY-X*r},p2:{x:A.pointX+L*r,y:A.pointY+X*r}},fitPoints:R.map(E=>({x:E.x,y:E.y}))}}function Ce(n){const t=Math.max(0,Math.min(1,n));return t<=.04045?t/12.92:Math.pow((t+.055)/1.055,2.4)}function ir(n,t,e){if(t<=0||e<=0||n.length!==t*e)return new Uint8Array(Math.max(0,t*e));let i=1/0,r=-1/0;for(let f=0;f<n.length;f++){const p=n[f];Number.isFinite(p)&&(p<i&&(i=p),p>r&&(r=p))}if(!Number.isFinite(i)||!Number.isFinite(r)||r<=i+1e-9)return new Uint8Array(t*e);const s=1024,o=new Uint32Array(s),a=r-i;for(let f=0;f<n.length;f++){const p=Math.max(0,Math.min(1,(n[f]-i)/a)),g=Math.min(s-1,Math.max(0,Math.floor(p*(s-1))));o[g]++}const l=n.length,u=f=>{const p=l*f;let g=0;for(let x=0;x<s;x++)if(g+=o[x],g>=p)return i+x/Math.max(1,s-1)*a;return r},d=u(.01),c=u(.99),m=Math.max(1e-9,c-d),h=new Uint8Array(t*e);for(let f=0;f<n.length;f++){const p=Math.max(0,Math.min(1,(n[f]-d)/m));h[f]=Math.round(p*255)}return h}function Hn(n,t,e=0){const i=new Float32Array(n.width*n.height),r=n.data;for(let s=0,o=0;s<r.length;s+=4,o++)i[o]=Fs(r,s,t,e);return ir(i,n.width,n.height)}function ul(n){return Number.isFinite(n)?Math.max(0,Math.min(65535,Number(n))):0}function Fs(n,t,e,i=0){let r=n[t]/255,s=n[t+1]/255,o=n[t+2]/255;e&&(r=Ce(r),s=Ce(s),o=Ce(o));const a=.2126*r+.7152*s+.0722*o;return Math.max(0,a-ul(i)/65535)}function hl(n,t){const e=n.width,i=n.height,r=n.data;if(r.length<e*i*3)return new Uint8Array(e*i);const s=new Float32Array(e*i);for(let o=0;o<e*i;o++){const a=o*3;t!==void 0?s[o]=r[a+t]:s[o]=.2126*r[a]+.7152*r[a+1]+.0722*r[a+2]}return ir(s,e,i)}function dl(n){const t=new Float32Array(n.width*n.height);for(let e=0;e<n.data.length;e++)t[e]=n.data[e];return ir(t,n.width,n.height)}function Yr(n,t,e){const i=Mt(t,n.width,n.height);if(!i)return null;const r=new Uint16Array(i.w*i.h*3),s=n.data;let o=0;for(let a=i.y;a<i.y+i.h;a++)for(let l=i.x;l<i.x+i.w;l++){const u=(a*n.width+l)*4;let d=s[u]/255,c=s[u+1]/255,m=s[u+2]/255;e&&(d=Ce(d),c=Ce(c),m=Ce(m)),r[o++]=Math.max(0,Math.min(65535,Math.round(d*65535))),r[o++]=Math.max(0,Math.min(65535,Math.round(c*65535))),r[o++]=Math.max(0,Math.min(65535,Math.round(m*65535)))}return{data:r,width:i.w,height:i.h}}function Wr(n,t,e,i=0){const r=Mt(t,n.width,n.height);if(!r)return null;const s=new Uint16Array(r.w*r.h*3),o=n.data;let a=0;for(let l=r.y;l<r.y+r.h;l++)for(let u=r.x;u<r.x+r.w;u++){const d=(l*n.width+u)*4,c=Math.max(0,Math.min(65535,Math.round(Fs(o,d,e,i)*65535)));s[a++]=c,s[a++]=c,s[a++]=c}return{data:s,width:r.w,height:r.h}}function fl(n,t){const e=Mt(t,n.width,n.height);if(!e)return null;const i=new Uint16Array(e.w*e.h),r=n.data;let s=0;for(let o=e.y;o<e.y+e.h;o++)for(let a=e.x;a<e.x+e.w;a++){const l=(o*n.width+a)*4;i[s++]=Math.max(0,Math.min(65535,Math.round((.2126*r[l]+.7152*r[l+1]+.0722*r[l+2])*257)))}return{data:i,width:e.w,height:e.h}}function pl(n,t){const e=Mt(t,n.width,n.height);if(!e)return null;const i=new Uint16Array(e.w*e.h);let r=0;for(let s=e.y;s<e.y+e.h;s++){const o=s*n.width;for(let a=e.x;a<e.x+e.w;a++)i[r++]=n.data[o+a]}return{data:i,width:e.w,height:e.h}}function qt(n,t,e){return{x:n.x*t,y:n.y*e}}function ml(n,t,e){return{p1:qt(n.p1,t,e),p2:qt(n.p2,t,e)}}function Yi(n,t){const e=t(n);return{x:Number.isFinite(e.x)?e.x:n.x,y:Number.isFinite(e.y)?e.y:n.y}}function gl(n,t){return n.map(e=>Yi(e,t))}function ge(n,t,e,i=0,r=0){if(!n||n.length<8)return;const s=n.map(o=>({x:o.x*t-i,y:o.y*e-r})).filter(o=>Number.isFinite(o.x)&&Number.isFinite(o.y));return s.length>=8?s:void 0}function yl(n,t){return{p1:Yi(n.p1,t),p2:Yi(n.p2,t)}}function cn(n,t,e){return{p1:{x:n.p1.x-t,y:n.p1.y-e},p2:{x:n.p2.x-t,y:n.p2.y-e}}}function xl(n,t,e,i){const r=Math.max(0,Math.min(n.width-1,t)),o=(Math.max(0,Math.min(n.height-1,e))*n.width+r)*4;let a=n.data[o]/255,l=n.data[o+1]/255,u=n.data[o+2]/255;return i&&(a=Ce(a),l=Ce(l),u=Ce(u)),(.2126*a+.7152*l+.0722*u)*65535}function As(n){return n.kind==="u16-mono"}function en(n){return n.width}function nn(n){return n.height}function ni(n,t,e,i){if(As(n)){const r=Math.max(0,Math.min(n.width-1,t)),s=Math.max(0,Math.min(n.height-1,e));return n.data[s*n.width+r]}return xl(n,t,e,i)}function bl(n,t,e,i){if(As(n)&&n.coordinateSpace==="distorted-padded"){const r=Math.round(n.paddingOffsetX??0),s=Math.round(n.paddingOffsetY??0);return ni(n,t+r,e+s,i)}return ni(n,t,e,i)}function _l(n,t,e,i,r=3){const o=[...n,{x:(n[0].x+n[1].x+n[2].x+n[3].x)*.25,y:(n[0].y+n[1].y+n[2].y+n[3].y)*.25},{x:(n[0].x+n[1].x)*.5,y:(n[0].y+n[1].y)*.5},{x:(n[1].x+n[2].x)*.5,y:(n[1].y+n[2].y)*.5},{x:(n[2].x+n[3].x)*.5,y:(n[2].y+n[3].y)*.5},{x:(n[3].x+n[0].x)*.5,y:(n[3].y+n[0].y)*.5}].map(a=>ke(a,t)).filter(a=>Number.isFinite(a.x)&&Number.isFinite(a.y));return o.length===0?null:Mt(Vt(o,r),e,i)}function rr(n,t,e,i){const r=new Map;for(let s=n.y;s<n.y+n.h;s++)for(let o=n.x;o<n.x+n.w;o++){const a=Ut({x:o,y:s},t);if(!Number.isFinite(a.x)||!Number.isFinite(a.y))continue;const l=Math.round(a.x),u=Math.round(a.y);if(l<0||u<0||l>=e||u>=i)continue;const d=r.get(u);d?(l<d.start&&(d.start=l),l>d.end&&(d.end=l)):r.set(u,{start:l,end:l})}return r}function Ml(n,t,e,i,r,s,o){const a=new Map,l=t.p2.x-t.p1.x,u=t.p2.y-t.p1.y,d=Math.hypot(l,u);if(!Number.isFinite(d)||d<=1e-6)return a;const c=l/d,m=u/d,h=-m,f=c,p={x:(t.p1.x+t.p2.x)*.5,y:(t.p1.y+t.p2.y)*.5},g=Math.max(1,e+1),x=Math.max(1,i+1.5);for(let y=n.y;y<n.y+n.h;y++)for(let b=n.x;b<n.x+n.w;b++){const _=b+.5,F=y+.5,v=_-p.x,w=F-p.y,P=v*c+w*m;if(!Number.isFinite(P)||Math.abs(P)>g)continue;const T=v*h+w*f;if(!Number.isFinite(T)||Math.abs(T)>x)continue;const C=Ut({x:_,y:F},r);if(!Number.isFinite(C.x)||!Number.isFinite(C.y))continue;const k=Math.round(C.x),M=Math.round(C.y);if(k<0||M<0||k>=s||M>=o)continue;const S=a.get(M);S?(k<S.start&&(S.start=k),k>S.end&&(S.end=k)):a.set(M,{start:k,end:k})}return a}function Ts(n,t,e,i,r,s){const o=new Map,a=t.p2.x-t.p1.x,l=t.p2.y-t.p1.y,u=Math.hypot(a,l);if(!Number.isFinite(u)||u<=1e-6)return o;const d=a/u,c=l/u,m=-c,h=d,f={x:(t.p1.x+t.p2.x)*.5,y:(t.p1.y+t.p2.y)*.5},p=Math.max(1,e+1),g=Math.max(1,i+1.5),x=Mt(n,r,s);if(!x)return o;for(let y=x.y;y<x.y+x.h;y++)for(let b=x.x;b<x.x+x.w;b++){const _=b-f.x,F=y-f.y,v=_*d+F*c;if(!Number.isFinite(v)||Math.abs(v)>p)continue;const w=_*m+F*h;if(!Number.isFinite(w)||Math.abs(w)>g)continue;const P=o.get(y);P?(b<P.start&&(P.start=b),b>P.end&&(P.end=b)):o.set(y,{start:b,end:b})}return o}function Is(n,t,e,i){const r=new Map;for(const[s,o]of n)for(let a=o.start;a<=o.end;a++){const l=ke({x:a,y:s},t);if(!Number.isFinite(l.x)||!Number.isFinite(l.y))continue;const u=Math.round(l.x),d=Math.round(l.y);if(u<0||d<0||u>=e||d>=i)continue;const c=r.get(d);c?(u<c.start&&(c.start=u),u>c.end&&(c.end=u)):r.set(d,{start:u,end:u})}return r}function sr(n){return Math.abs(n.k1)<1e-4&&Math.abs(n.k2)<1e-4}function wl(n){return[{x:n.x,y:n.y},{x:n.x+n.w,y:n.y},{x:n.x+n.w,y:n.y+n.h},{x:n.x,y:n.y+n.h}]}function Kt(n,t,e,i,r){return ke({x:i.x+n*t,y:i.y+n*e},r)}function ar(n,t,e,i,r){const o=Kt(n,t,e,i,r),a=Kt(n+1e-4,t,e,i,r);return{x:(a.x-o.x)/1e-4,y:(a.y-o.y)/1e-4}}function Rs(n,t,e,i,r,s){let o=.01;const a=d=>{const c=Kt(d,t,e,i,s);return Math.hypot(c.x-r.x,c.y-r.y)},l=a(n),u=a(n+o);if(!Number.isFinite(l)||!Number.isFinite(u))return null;if(l>u){let d=n,c=n+o;for(let m=0;m<24;m++){o*=2;const h=d+o,f=a(h),p=a(c);if(!Number.isFinite(f)||!Number.isFinite(p))break;if(f>=p)return{a:d,b:h};d=c,c=h}}else{let d=n,c=n+o;for(let m=0;m<24;m++){o*=2;const h=c-o,f=a(h),p=a(d);if(!Number.isFinite(f)||!Number.isFinite(p))break;if(f>=p)return{a:h,b:c};c=d,d=h}}return{a:n-Math.max(.5,o),b:n+Math.max(.5,o)}}function Sl(n,t,e=33){const i=n.p2.x-n.p1.x,r=n.p2.y-n.p1.y,s=Math.hypot(i,r);if(!Number.isFinite(s)||s<=1e-6)return[ke(n.p1,t),ke(n.p2,t)];const o=i/s,a=r/s,l={x:(n.p1.x+n.p2.x)*.5,y:(n.p1.y+n.p2.y)*.5},u=s*.5,d=Math.max(9,e),c=[];for(let m=0;m<d;m++){const h=d===1?.5:m/(d-1),f=-u+h*(u*2);c.push(Kt(f,o,a,l,t))}return c}function vl(n,t,e,i,r,s,o=1){const a=n.p2.x-n.p1.x,l=n.p2.y-n.p1.y,u=Math.hypot(a,l);if(!Number.isFinite(u)||u<=1e-6)return null;const d=a/u,c=l/u,m={x:(n.p1.x+n.p2.x)*.5,y:(n.p1.y+n.p2.y)*.5},h=Math.max(24,Math.round(e*2)+1),f=[];for(let p=0;p<h;p++){const g=h===1?.5:p/(h-1),x=-e+g*(e*2),y=Kt(x,d,c,m,t),b=ar(x,d,c,m,t),_=Math.hypot(b.x,b.y);if(!Number.isFinite(_)||_<=1e-9)continue;const F=-b.y/_,v=b.x/_;f.push({x:y.x+F*(i+o),y:y.y+v*(i+o)},{x:y.x-F*(i+o),y:y.y-v*(i+o)})}if(f.length<2){const p={p1:ke(n.p1,t),p2:ke(n.p2,t)},g=be(p,i+o);return g?Mt(Vt(g,o),r,s):null}return Mt(Vt(f,o),r,s)}function tn(n,t,e,i,r){return Ut({x:i.x+n*t,y:i.y+n*e},r)}function Pl(n,t,e,i,r,s){let o=.01;const a=d=>{const c=tn(d,t,e,i,s);return Math.hypot(c.x-r.x,c.y-r.y)},l=a(n),u=a(n+o);if(!Number.isFinite(l)||!Number.isFinite(u))return null;if(l>u){let d=n,c=n+o;for(let m=0;m<24;m++){o*=2;const h=d+o,f=a(h),p=a(c);if(!Number.isFinite(f)||!Number.isFinite(p))break;if(f>=p)return{a:d,b:h};d=c,c=h}}else{let d=n,c=n+o;for(let m=0;m<24;m++){o*=2;const h=c-o,f=a(h),p=a(d);if(!Number.isFinite(f)||!Number.isFinite(p))break;if(f>=p)return{a:h,b:c};c=d,d=h}}return{a:n-Math.max(.5,o),b:n+Math.max(.5,o)}}function or(n,t,e){const i=(t.x-n.x)*(t.y-e.y)-(t.x-e.x)*(t.y-n.y);if(Math.abs(i)<=1e-12)return .5*(n.x+e.x);const r=(t.x-n.x)*(t.x-n.x)*(t.y-e.y)-(t.x-e.x)*(t.x-e.x)*(t.y-n.y),s=t.x-.5*r/i;return Number.isFinite(s)?s:.5*(n.x+e.x)}function Cl(n,t,e,i,r){const o=tn(n,t,e,i,r),a=tn(n+1e-4,t,e,i,r);return{x:(a.x-o.x)/1e-4,y:(a.y-o.y)/1e-4}}function kl(n,t,e,i,r,s,o=!1,a){const l=[],u=[],d=a?Lt*2:Lt,c=Math.max(1,Math.min(s,d)),m=i.p2.x-i.p1.x,h=i.p2.y-i.p1.y,f=Math.hypot(m,h);if(!Number.isFinite(f)||f<=1e-6)return null;const p=m/f,g=h/f,x={x:(i.p1.x+i.p2.x)*.5,y:(i.p1.y+i.p2.y)*.5},y={p1:Ut(i.p1,e),p2:Ut(i.p2,e)},b=y.p2.x-y.p1.x,_=y.p2.y-y.p1.y,F=Math.hypot(b,_);if(!Number.isFinite(F)||F<=1e-6)return null;const v=b/F,w=_/F,P=-w,T=v,C={x:(y.p1.x+y.p2.x)*.5,y:(y.p1.y+y.p2.y)*.5},k=Mt(a||Vt(be(i,s+2)??[i.p1,i.p2],2),n.width,n.height);if(!k)return null;const M=Ml(k,i,r,c,e,en(t),nn(t));if(M.size===0)return null;const S=!sr(e);for(const[E,U]of M)if(!(E<0||E>=nn(t)))for(let O=U.start;O<=U.end;O++){if(O<0||O>=en(t))continue;const z={x:O,y:E};let J,V;if(S){const Q=ke(z,e);if(!Number.isFinite(Q.x)||!Number.isFinite(Q.y)||Math.round(Q.x)<0||Math.round(Q.x)>=n.width||Math.round(Q.y)<0||Math.round(Q.y)>=n.height)continue;const W=Q.x-x.x,tt=Q.y-x.y,at=W*p+tt*g;if(!Number.isFinite(at))continue;J=at,V=W*-g+tt*p;const B=Pl(at,p,g,x,z,e);if(!B)continue;const Z=.5*(B.a+B.b),q=tn(B.a,p,g,x,e),I=tn(Z,p,g,x,e),H=tn(B.b,p,g,x,e),G=or({x:B.a,y:Math.hypot(q.x-z.x,q.y-z.y)},{x:Z,y:Math.hypot(I.x-z.x,I.y-z.y)},{x:B.b,y:Math.hypot(H.x-z.x,H.y-z.y)});if(!Number.isFinite(G))continue;J=G;const Y=Cl(G,p,g,x,e),et=Math.hypot(Y.x,Y.y);if(!Number.isFinite(et)||et<=1e-9)continue;const nt=Y.x/et,ct=-(Y.y/et),dt=nt,rt=tn(G,p,g,x,e);V=(z.x-rt.x)*ct+(z.y-rt.y)*dt}else{const Q=z.x-C.x,W=z.y-C.y;J=Q*v+W*w,V=Q*P+W*T}!Number.isFinite(J)||Math.abs(J)>r||!Number.isFinite(V)||Math.abs(V)>c||(l.push(V),u.push(ni(t,O,E,o)))}if(l.length<8)return null;const R=Ut(i.p1,e),A=Ut(i.p2,e),L=A.x-R.x,X=A.y-R.y,D=Math.abs(L)>=Math.abs(X)?1:2;return xn(l,u,D,d)}function Fl(n,t,e,i,r,s,o=!1,a,l,u,d=!1){const c=[],m=[],h=a?Lt*2:Lt,f=Math.max(1,Math.min(s,h)),p=i.p2.x-i.p1.x,g=i.p2.y-i.p1.y,x=Math.hypot(p,g);if(!Number.isFinite(x)||x<=1e-6)return null;const y=p/x,b=g/x,_=-b,F=y,v={x:(i.p1.x+i.p2.x)*.5,y:(i.p1.y+i.p2.y)*.5},w=Mt(a||Vt(be(i,h*4+2)??[i.p1,i.p2],2),n.width,n.height);if(!w)return null;const P=l??(u?rr(Mt(u,en(t),nn(t))??u,e,n.width,n.height):Ts(w,i,Math.max(1,r),f*4+.5,n.width,n.height));if(P.size===0)return null;const T=Is(P,e,en(t),nn(t));if(T.size===0)return null;const C=!sr(e);for(const[M,S]of T)for(let R=S.start;R<=S.end;R++){const A={x:R,y:M},L=Ut(A,e);if(!Number.isFinite(L.x)||!Number.isFinite(L.y)||Math.round(L.x)<0||Math.round(L.x)>=n.width||Math.round(L.y)<0||Math.round(L.y)>=n.height)continue;const X=L.x-v.x,D=L.y-v.y,E=X*y+D*b;let U=X*_+D*F;if(C){const O=Rs(E,y,b,v,A,e);if(!O)continue;const z=.5*(O.a+O.b),J=Kt(O.a,y,b,v,e),V=Kt(z,y,b,v,e),Q=Kt(O.b,y,b,v,e),W=or({x:O.a,y:Math.hypot(J.x-A.x,J.y-A.y)},{x:z,y:Math.hypot(V.x-A.x,V.y-A.y)},{x:O.b,y:Math.hypot(Q.x-A.x,Q.y-A.y)});if(!Number.isFinite(W))continue;const tt=ar(W,y,b,v,e),at=Math.hypot(tt.x,tt.y);if(!Number.isFinite(at)||at<=1e-9)continue;const B=tt.x/at,q=-(tt.y/at),I=B,H=Kt(W,y,b,v,e);U=(A.x-H.x)*q+(A.y-H.y)*I}!Number.isFinite(E)||Math.abs(E)>Math.max(1,r)||!Number.isFinite(U)||Math.abs(U)>f||(c.push(U),m.push(bl(t,R,M,o)))}if(c.length<8)return null;const k=Math.abs(p)>=Math.abs(g)?1:2;return d?Fn(c,m,k,h):xn(c,m,k,h)}function Al(n,t,e,i,r,s){const o=n.width,a=n.height,l=(n.bayerPattern||"RGGB").toUpperCase(),u=s!=null&&s.correctedRect?Lt*2:Lt,d=Math.max(1,Math.min(r,u)),c=(s==null?void 0:s.restrictToStrip)??!0,m=e.p2.x-e.p1.x,h=e.p2.y-e.p1.y,f=Math.hypot(m,h);if(!Number.isFinite(f)||f<=1e-6)return null;const p=m/f,g=h/f,x=-g,y=p,b={x:(e.p1.x+e.p2.x)*.5,y:(e.p1.y+e.p2.y)*.5},_={p1:{x:b.x-p*Math.max(1,i),y:b.y-g*Math.max(1,i)},p2:{x:b.x+p*Math.max(1,i),y:b.y+g*Math.max(1,i)}},F=be(_,d+2),v=(s!=null&&s.fixedRawRect?Mt(s.fixedRawRect,o,a):null)??(s!=null&&s.correctedRect?er(s.correctedRect,t,o,a):null)??(F?_l(F,t,o,a,2):null);if(!v)return null;const w=[],P=[];for(let k=v.y;k<v.y+v.h;k++){const M=k*o;for(let S=v.x;S<v.x+v.w;S++){if(!vt(S,k,l,s==null?void 0:s.greenPhase))continue;const R=Ut({x:S,y:k},t);if(!Number.isFinite(R.x)||!Number.isFinite(R.y))continue;const A=R.x-b.x,L=R.y-b.y,X=A*p+L*g;if(!Number.isFinite(X)||c&&Math.abs(X)>Math.max(1,i))continue;const D=A*x+L*y;if(!Number.isFinite(D)||c&&Math.abs(D)>d)continue;w.push(D);let E;E=Math.max(0,n.data[M+S]-yn(s==null?void 0:s.blackLevel,S,k)),P.push(E)}}if(w.length<8)return null;const T=Math.abs(m)>=Math.abs(h)?1:2,C=Math.max(2,(s==null?void 0:s.shortSidePxOverride)??(c?d*2:Math.min(v.w,v.h)));return ui(w,P,C,s==null?void 0:s.manualBinSize,T,s==null?void 0:s.preferAutoPerEdgeBin,!1,!!(s!=null&&s.forceLegacyModel))}function Ns(n,t,e,i,r,s=!1,o,a,l,u=!1){if(a&&l){const P=l.p2.x-l.p1.x,T=l.p2.y-l.p1.y,C=Math.hypot(P,T);if(Number.isFinite(C)&&C>1e-6)return kl(a,n,t,l,Math.max(1,C*.5),r,s,o)}const d=[],c=[],m=Lt,h=e.p2.x-e.p1.x,f=e.p2.y-e.p1.y,p=Math.hypot(h,f);if(!Number.isFinite(p)||p<=1e-6)return null;const g=h/p,x=f/p,y={x:(e.p1.x+e.p2.x)*.5,y:(e.p1.y+e.p2.y)*.5},b=(o?Mt(o,en(n),nn(n)):null)??vl(e,t,i+1,r+1,en(n),nn(n),1);if(!b)return null;const _=rr(b,t,en(n),nn(n));if(_.size===0)return null;const F=-x,v=g;for(const[P,T]of _)for(let C=T.start;C<=T.end;C++){const k={x:C,y:P};let M=(k.x-y.x)*g+(k.y-y.y)*x,S=(k.x-y.x)*F+(k.y-y.y)*v;!Number.isFinite(M)||Math.abs(M)>i+1||!Number.isFinite(S)||Math.abs(S)>=m||(d.push(S),c.push(ni(n,C,P,s)))}if(d.length<8)return null;const w=Math.abs(h)>=Math.abs(f)?1:2;return u?Fn(d,c,w,m):xn(d,c,w,m)}function Tl(n,t,e){var m;const i=e.sourceMode??(t.isThreePlane?"three-plane":"rggb-raw"),r=e.useQuadraticProjection!==!1,s=!!e.forceRenderedMeasurement,o=n.width,a=n.height,l=e.threePlaneChannel,u=nr(e.detectionTuning),d=e.monochromeBlackLevel??0;if(i==="rggb-raw"&&!s){if(!t||t.isThreePlane)return null;const h=Bo(t,e.greenPhase),f=o/Math.max(1,t.width),p=a/Math.max(1,t.height);return{sourceMode:i,detectionGray:h,detectionWidth:t.width,detectionHeight:t.height,detectToMeasureX:1,detectToMeasureY:1,detectToDisplayX:f,detectToDisplayY:p,measureToDisplayX:f,measureToDisplayY:p,detectPointToDisplay:g=>qt(g,f,p),measurePointToDisplay:g=>qt(g,f,p),displayPointToDetect:g=>qt(g,1/Math.max(1e-9,f),1/Math.max(1e-9,p)),measureUsesDisplayLine:!1,measureWidth:t.width,measureHeight:t.height,refineLine:(g,x,y)=>cl(t.data,t.width,t.height,{p1:g,p2:x},y*.5,Math.max(4,y*.2),{greenPhase:e.greenPhase,bayerPattern:t.bayerPattern})||ye(h,t.width,t.height,g,x,y*.5,Math.max(4,y*.2)),measureEdge:(g,x,y,b,_)=>ve(t.data,t.width,t.height,g,x,y,b,{greenOnly:!0,greenPhase:e.greenPhase,bayerPattern:t.bayerPattern,blackLevel:e.blackLevel??void 0,preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?ge(_==null?void 0:_.fitPoints,1,1):void 0})}}if(s){const h=!!e.distortionCurveApplied&&!!e.distortionModel,f=i==="rggb-raw"&&!!e.distortionCorrected&&!!e.distortionModel&&!t.isThreePlane,p=!!e.distortionCorrected&&!!e.distortionModel&&!!e.distortionOriginalSamplingPlane,g=!!e.distortionCorrected&&!!e.distortionSamplingPlane,x=n,y=Hn(x,!!e.sfrHasGamma,i==="unmix-bw"?d:0);return{sourceMode:i,detectionGray:y,detectionWidth:x.width,detectionHeight:x.height,detectToMeasureX:1,detectToMeasureY:1,detectToDisplayX:1,detectToDisplayY:1,measureToDisplayX:1,measureToDisplayY:1,detectPointToDisplay:b=>b,measurePointToDisplay:b=>b,displayPointToDetect:b=>b,measureUsesDisplayLine:!1,measureWidth:x.width,measureHeight:x.height,refineLine:(b,_,F)=>(p?ol(y,n.width,n.height,b,_,Lt):null)||ye(y,n.width,n.height,b,_,F*.5,Math.max(4,F*.2)),measureEdge:(b,_,F,v,w)=>{const P=e.distortionModel?er(b,e.distortionModel,t.width,t.height):null;if(f){const M={p1:Ut(_.p1,e.distortionModel),p2:Ut(_.p2,e.distortionModel)},S=Math.hypot(M.p2.x-M.p1.x,M.p2.y-M.p1.y),R=Math.max(2,S*.5*u.sampleHalfWidthRatio);return eo(t,e.distortionModel,M,Math.max(1,S*.5),R,{greenPhase:e.greenPhase,blackLevel:e.blackLevel??void 0,correctedRect:b})}if(f)return Al(t,e.distortionModel,_,F,v,{greenPhase:e.greenPhase,blackLevel:e.blackLevel??void 0,correctedRect:b,fixedRawRect:P,preferAutoPerEdgeBin:!0});if(p)return Fl(n,e.distortionOriginalSamplingPlane,e.distortionModel,_,F,v,!!e.sfrHasGamma,b,(w==null?void 0:w.correctedScanlines)??null,P);if(h){const M={p1:Ut(_.p1,e.distortionModel),p2:Ut(_.p2,e.distortionModel)},S=Math.hypot(M.p2.x-M.p1.x,M.p2.y-M.p1.y);return Ns(e.distortionOriginalSamplingPlane??e.distortionSamplingPlane??e.distortionSamplingImage??n,e.distortionModel,M,Math.max(1,S*.5),v,!!e.sfrHasGamma,b,e.distortionBaseImage??n,_)}if(g){const M=pl(e.distortionSamplingPlane,b);if(!M)return null;const S=cn(_,b.x,b.y);return ve(M.data,M.width,M.height,{x:0,y:0,w:M.width,h:M.height},S,F,v,{preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?ge(w==null?void 0:w.fitPoints,1,1,b.x,b.y):void 0})}const C=i==="unmix-bw"?Wr(n,b,!!e.sfrHasGamma,d):Yr(n,b,!!e.sfrHasGamma);if(!C)return null;const k=cn(_,b.x,b.y);return ve(C.data,C.width,C.height,{x:0,y:0,w:C.width,h:C.height},k,F,v,{isThreePlane:!0,threePlaneChannel:void 0,preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?ge(w==null?void 0:w.fitPoints,1,1,b.x,b.y):void 0})}}}if(i==="three-plane"){if(t.isThreePlane&&!e.sfrHasGamma){const f=hl(t,l),p=o/Math.max(1,t.width),g=a/Math.max(1,t.height);return{sourceMode:i,detectionGray:f,detectionWidth:t.width,detectionHeight:t.height,detectToMeasureX:1,detectToMeasureY:1,detectToDisplayX:p,detectToDisplayY:g,measureToDisplayX:p,measureToDisplayY:g,detectPointToDisplay:x=>qt(x,p,g),measurePointToDisplay:x=>qt(x,p,g),displayPointToDetect:x=>qt(x,1/Math.max(1e-9,p),1/Math.max(1e-9,g)),measureUsesDisplayLine:!1,measureWidth:t.width,measureHeight:t.height,refineLine:(x,y,b)=>ye(f,t.width,t.height,x,y,b*.5,Math.max(4,b*.2)),measureEdge:(x,y,b,_,F)=>ve(t.data,t.width,t.height,x,y,b,_,{isThreePlane:!0,threePlaneChannel:l,preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?ge(F==null?void 0:F.fitPoints,1,1):void 0})}}const h=Hn(n,!!e.sfrHasGamma);return{sourceMode:i,detectionGray:h,detectionWidth:n.width,detectionHeight:n.height,detectToMeasureX:1,detectToMeasureY:1,detectToDisplayX:1,detectToDisplayY:1,measureToDisplayX:1,measureToDisplayY:1,detectPointToDisplay:f=>f,measurePointToDisplay:f=>f,displayPointToDetect:f=>f,measureUsesDisplayLine:!1,measureWidth:n.width,measureHeight:n.height,refineLine:(f,p,g)=>ye(h,n.width,n.height,f,p,g*.5,Math.max(4,g*.2)),measureEdge:(f,p,g,x,y)=>{const b=Yr(n,f,!!e.sfrHasGamma);if(!b)return null;const _=cn(p,f.x,f.y);return ve(b.data,b.width,b.height,{x:0,y:0,w:b.width,h:b.height},_,g,x,{isThreePlane:!0,threePlaneChannel:l,preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?ge(y==null?void 0:y.fitPoints,1,1,f.x,f.y):void 0})}}}if(i==="unmix-bw"){if(t&&!t.isThreePlane&&e.displaySettings){const f=Da(t,e.displaySettings,e.blackLevel??e.monochromeBlackLevel??void 0);if(f){const p=dl(f),g=o/Math.max(1,t.width),x=a/Math.max(1,t.height);return{sourceMode:i,detectionGray:p,detectionWidth:t.width,detectionHeight:t.height,detectToMeasureX:1,detectToMeasureY:1,detectToDisplayX:g,detectToDisplayY:x,measureToDisplayX:g,measureToDisplayY:x,detectPointToDisplay:y=>qt(y,g,x),measurePointToDisplay:y=>qt(y,g,x),displayPointToDetect:y=>qt(y,1/Math.max(1e-9,g),1/Math.max(1e-9,x)),measureUsesDisplayLine:!1,measureWidth:t.width,measureHeight:t.height,refineLine:(y,b,_)=>ye(p,t.width,t.height,y,b,_*.5,Math.max(4,_*.2)),measureEdge:(y,b,_,F,v)=>ve(f.data,f.width,f.height,y,b,_,F,{preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?ge(v==null?void 0:v.fitPoints,1,1):void 0})}}}const h=Hn(n,!!e.sfrHasGamma,d);return{sourceMode:i,detectionGray:h,detectionWidth:n.width,detectionHeight:n.height,detectToMeasureX:1,detectToMeasureY:1,detectToDisplayX:1,detectToDisplayY:1,measureToDisplayX:1,measureToDisplayY:1,detectPointToDisplay:f=>f,measurePointToDisplay:f=>f,displayPointToDetect:f=>f,measureUsesDisplayLine:!1,measureWidth:n.width,measureHeight:n.height,refineLine:(f,p,g)=>ye(h,n.width,n.height,f,p,g*.5,Math.max(4,g*.2)),measureEdge:(f,p,g,x,y)=>{const b=Wr(n,f,!!e.sfrHasGamma,d);if(!b)return null;const _=cn(p,f.x,f.y);return ve(b.data,b.width,b.height,{x:0,y:0,w:b.width,h:b.height},_,g,x,{isThreePlane:!0,threePlaneChannel:void 0,preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?ge(y==null?void 0:y.fitPoints,1,1,f.x,f.y):void 0})}}}const c=Hn(n,!1);if(t&&!t.isThreePlane&&((m=e.displaySettings)==null?void 0:m.renderMode)==="advanced-zero-dep"&&e.displaySettings.advancedZeroDep){const h=o/Math.max(1,t.width),f=a/Math.max(1,t.height);return{sourceMode:i,detectionGray:c,detectionWidth:n.width,detectionHeight:n.height,detectToMeasureX:t.width/Math.max(1,n.width),detectToMeasureY:t.height/Math.max(1,n.height),detectToDisplayX:1,detectToDisplayY:1,measureToDisplayX:h,measureToDisplayY:f,detectPointToDisplay:p=>p,measurePointToDisplay:p=>qt(p,h,f),displayPointToDetect:p=>p,measureUsesDisplayLine:!1,measureWidth:t.width,measureHeight:t.height,refineLine:(p,g,x)=>ye(c,n.width,n.height,p,g,x*.5,Math.max(4,x*.2)),measureEdge:(p,g,x,y,b)=>{const _=Ua(t,p,e.displaySettings);if(!_||_.width<8||_.height<8)return null;const F=cn(g,p.x,p.y);return ve(_.data,_.width,_.height,{x:0,y:0,w:_.width,h:_.height},F,x,y,{preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?ge(b==null?void 0:b.fitPoints,t.width/Math.max(1,n.width),t.height/Math.max(1,n.height),p.x,p.y):void 0})}}}if(t&&!t.isThreePlane){const h=o/Math.max(1,t.width),f=a/Math.max(1,t.height);return{sourceMode:i,detectionGray:c,detectionWidth:n.width,detectionHeight:n.height,detectToMeasureX:t.width/Math.max(1,n.width),detectToMeasureY:t.height/Math.max(1,n.height),detectToDisplayX:1,detectToDisplayY:1,measureToDisplayX:h,measureToDisplayY:f,detectPointToDisplay:p=>p,measurePointToDisplay:p=>qt(p,h,f),displayPointToDetect:p=>p,measureUsesDisplayLine:!1,measureWidth:t.width,measureHeight:t.height,refineLine:(p,g,x)=>ye(c,n.width,n.height,p,g,x*.5,Math.max(4,x*.2)),measureEdge:(p,g,x,y,b)=>ve(t.data,t.width,t.height,p,g,x,y,{blackLevel:e.blackLevel??void 0,preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?ge(b==null?void 0:b.fitPoints,t.width/Math.max(1,n.width),t.height/Math.max(1,n.height)):void 0})}}return{sourceMode:i,detectionGray:c,detectionWidth:n.width,detectionHeight:n.height,detectToMeasureX:1,detectToMeasureY:1,detectToDisplayX:1,detectToDisplayY:1,measureToDisplayX:1,measureToDisplayY:1,detectPointToDisplay:h=>h,measurePointToDisplay:h=>h,displayPointToDetect:h=>h,measureUsesDisplayLine:!1,measureWidth:n.width,measureHeight:n.height,refineLine:(h,f,p)=>ye(c,n.width,n.height,h,f,p*.5,Math.max(4,p*.2)),measureEdge:(h,f,p,g,x)=>{const y=fl(n,h);if(!y)return null;const b=cn(f,h.x,h.y);return To(y.data,y.width,y.height,{x:0,y:0,w:y.width,h:y.height},b,p,g,{blackLevel:e.monochromeBlackLevel??void 0,preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?ge(x==null?void 0:x.fitPoints,1,1,h.x,h.y):void 0})}}}function Il(n,t,e){var u,d,c,m,h,f,p,g;if(!n||!t)return[];(u=e.onProgress)==null||u.call(e,"Preparing source context...",0);const i=Tl(n,t,e);if(!i)return[];(d=e.onProgress)==null||d.call(e,"Preparing source context...",.08);const r=nr(e.detectionTuning),s=Math.min(1e3,Math.max(1,e.maxRegions??1e3)),o=Math.max(4,e.maxEdges??s*4);(c=e.onProgress)==null||c.call(e,"Detecting candidates...",.12);const a=Jo(i.detectionGray,i.detectionWidth,i.detectionHeight,s,e.detectionTuning,(x,y)=>{var b;(b=e.onProgress)==null||b.call(e,x,.12+.08*Math.max(0,Math.min(1,y)))});(m=e.onProgress)==null||m.call(e,"Detecting candidates...",.2);const l=[];for(let x=0;x<a.length;x++){const y=a.length<=0?1:x/a.length;if((h=e.onProgress)==null||h.call(e,`Measuring edges: region ${x+1}/${a.length}`,.2+.72*Math.min(1,y)),l.length>=o)break;const b=a[x],_=b.corners,F=`auto-region-${x+1}`;for(let v=0;v<4&&((f=e.onProgress)==null||f.call(e,`Measuring edges: region ${x+1}/${a.length}, edge ${v+1}/4`,.2+.72*Math.min(1,(x+v/4)/Math.max(1,a.length))),!(l.length>=o));v=v+1){const w=_[v],P=_[(v+1)%4],T=P.x-w.x,C=P.y-w.y,k=Math.hypot(T,C);if(!Number.isFinite(k)||k<24)continue;const M=.125,S={x:w.x+T*M,y:w.y+C*M},R={x:P.x-T*M,y:P.y-C*M},A=Math.hypot(R.x-S.x,R.y-S.y);if(!Number.isFinite(A)||A<12)continue;const L=i.refineLine(S,R,A),X=(L!=null&&L.fitPoints?Ur(L.fitPoints):null)||(L==null?void 0:L.line)||{p1:S,p2:R},D=ml(X,i.detectToMeasureX,i.detectToMeasureY),E=gl((L==null?void 0:L.fitPoints)??[],i.detectPointToDisplay),U=(E.length>=2?Ur(E):null)||yl(D,i.measurePointToDisplay),O=i.measureUsesDisplayLine?U:D,z=O.p2.x-O.p1.x,J=O.p2.y-O.p1.y,V=Math.hypot(z,J);if(!Number.isFinite(V)||V<=1e-6)continue;const Q=U.p2.x-U.p1.x,W=U.p2.y-U.p1.y,tt=Math.hypot(Q,W);if(!Number.isFinite(tt)||tt<=1e-6)continue;const at=!!e.distortionCurveApplied&&!!e.distortionModel,B=Q/tt,Z=W/tt;let q=Z,I=-B;const H=(U.p1.x+U.p2.x)*.5,G=(U.p1.y+U.p2.y)*.5,Y=i.detectPointToDisplay({x:b.centerX,y:b.centerY}),et=Y.x,nt=Y.y;(H-et)*q+(G-nt)*I<0&&(q=-q,I=-I);const it=V*.5,ct=Math.max(2,V*r.sampleHalfWidthRatio),dt=Math.max(2,tt*r.sampleHalfWidthRatio),rt=at?{p1:Ut(U.p1,e.distortionModel),p2:Ut(U.p2,e.distortionModel)}:void 0,Dt=rt?Math.max(1,Math.hypot(rt.p2.x-rt.p1.x,rt.p2.y-rt.p1.y)*.5):it,Gt=at?U:D,j=at?dt:ct,ut=be(Gt,j);if(!ut)continue;const ot=Mt(Vt(ut,2),at?n.width:i.measureWidth,at?n.height:i.measureHeight);if(!ot)continue;const pt=e.distortionCorrected&&e.distortionModel&&i.sourceMode==="rggb-raw"?er(ot,e.distortionModel,t.width,t.height):null,mt=at?Ns(e.distortionSamplingPlane??e.distortionSamplingImage??n,e.distortionModel,rt,Dt,j,!!e.sfrHasGamma,ot,e.distortionBaseImage??n,Gt):i.measureEdge(ot,O,it,j,L);if(!mt||(mt.autoLikeUsed=!0,!Fo(mt,e.useDeshading,0)))continue;const ee=e.useNR?-1:12,Bt=Eo([mt],ee,null,e.useDeshading,0,!0);if(!Bt||Bt.mtf50===null||!Ao(Bt.lsfCropped))continue;const _t=rt?Sl(rt,e.distortionModel,Math.max(21,Math.round(tt*.5))):mt.quadraticProjectionUsed?ko(E,U,Math.max(21,Math.round(tt*.5))):void 0,ce=rt&&_t&&_t.length>=2?Vt(_t,dt+2):null,$t=ce?wl(ce):be(U,dt);if(!$t)continue;const Fe=e.distortionCorrected?ot:ce??Vt($t,2);let $={x:H+q*(dt+12),y:G+I*(dt+12)},Ft=Dr(B,Z);if(_t&&_t.length>=3){const Tt=Math.floor(_t.length/2),xt=_t[Math.max(0,Tt-1)],wt=_t[Math.min(_t.length-1,Tt+1)],ne=_t[Tt],Jt=wt.x-xt.x,bn=wt.y-xt.y,Be=Math.hypot(Jt,bn);if(Be>1e-6){const Oe=bn/Be,Ve=-Jt/Be;Ft=Dr(Jt/Be,bn/Be);const Ge={x:ne.x-et,y:ne.y-nt},Xe=Ge.x*Oe+Ge.y*Ve>=0?1:-1;$={x:ne.x+Oe*Xe*(dt+12),y:ne.y+Ve*Xe*(dt+12)}}}l.push({id:`${F}-edge-${v+1}`,regionId:F,sourceMode:i.sourceMode,edgeIndex:v,label:Bt.mtf50.toFixed(3),mtf50:Bt.mtf50,angle:Ft,orientation:mt.orientation,edgeData:mt,sourceRect:Fe,rawSourceRect:(i.sourceMode==="rggb-raw"?pt??ot:pt)??void 0,quad:$t,line:U,originalLine:D,curveBaseLine:rt,curvePoints:_t,labelPoint:$,ridgePoints:E,outerSideMeans:b.outerSideMeans,outerSideQuads:b.outerSideQuads,distortionCorrected:e.distortionCorrected??!1})}}return(p=e.onProgress)==null||p.call(e,"Finalizing results...",.98),(g=e.onProgress)==null||g.call(e,"Finalizing results...",1),l}const Rl=n=>!n.blackLevels||n.blackLevels.length<4?null:[Number(n.blackLevels[0])||0,Number(n.blackLevels[1])||0,Number(n.blackLevels[2])||0,Number(n.blackLevels[3])||0],Ui=(n,t)=>{t instanceof ArrayBuffer&&(n.includes(t)||n.push(t))};self.onmessage=async n=>{var o,a;const{id:t,buffer:e,detect:i,options:r}=n.data,s=performance.now();try{const l=performance.now(),u=await Ea(e),d=performance.now()-l;let c=0,m=[];if(i&&!u.isXTrans){const f=u.isThreePlane?"three-plane":"rggb-raw",p=performance.now();m=Il({width:u.width,height:u.height},u,{...r,sourceMode:f,forceRenderedMeasurement:!1,blackLevel:(r==null?void 0:r.blackLevel)??Rl(u),onProgress:(g,x)=>{self.postMessage({id:t,type:"progress",stage:g,progress:x})}}),c=performance.now()-p}const h=[];Ui(h,e),Ui(h,(o=u.data)==null?void 0:o.buffer),Ui(h,(a=u.floatData)==null?void 0:a.buffer),self.postMessage({id:t,type:"result",success:!0,raw:u,rawFileBuffer:e,measurements:m,timings:{decodeMs:d,detectMs:c,totalMs:performance.now()-s}},h)}catch(l){self.postMessage({id:t,type:"result",success:!1,error:(l==null?void 0:l.message)||String(l)})}};
