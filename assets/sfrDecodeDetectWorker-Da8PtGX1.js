var Ws=Object.defineProperty;var js=(n,t,e)=>t in n?Ws(n,t,{enumerable:!0,configurable:!0,writable:!0,value:e}):n[t]=e;var bt=(n,t,e)=>js(n,typeof t!="symbol"?t+"":t,e);function Hs(n){const[[t,e,i],[r,s,o],[a,l,u]]=n,h=t*(s*u-o*l)-e*(r*u-o*a)+i*(r*l-s*a);if(Math.abs(h)<1e-10)return null;const c=1/h;return[[(s*u-o*l)*c,(i*l-e*u)*c,(e*o-i*s)*c],[(o*a-r*u)*c,(t*u-i*a)*c,(i*r-t*o)*c],[(r*l-s*a)*c,(e*a-t*l)*c,(t*s-e*r)*c]]}function Qs(n,t){const[[e,i,r],[s,o,a],[l,u,h]]=n,[c,g,d]=t;return[e*c+i*g+r*d,s*c+o*g+a*d,l*c+u*g+h*d]}class In{constructor(t){bt(this,"size");bt(this,"isPowerOfTwo");bt(this,"_real");bt(this,"_imag");bt(this,"_scratch",null);bt(this,"_rev",null);bt(this,"_m",0);bt(this,"_internalFFT",null);bt(this,"_chirpReal",null);bt(this,"_chirpImag",null);bt(this,"_bReal",null);bt(this,"_bImag",null);bt(this,"_hanning",null);bt(this,"_windowSumSq",0);this.size=t,this.isPowerOfTwo=(t&t-1)===0&&t>0,this._real=new Float32Array(t),this._imag=new Float32Array(t),this.isPowerOfTwo?this.initRadix2():this.initBluestein()}initRadix2(){const t=this.size,e=Math.log2(t);this._rev=new Uint32Array(t);for(let i=0;i<t;i++){let r=0,s=i;for(let o=0;o<e;o++)r=r<<1|s&1,s>>>=1;this._rev[i]=r}}initBluestein(){const t=this.size;this._m=Math.pow(2,Math.ceil(Math.log2(2*t-1))),this._internalFFT=new In(this._m),this._chirpReal=new Float32Array(t),this._chirpImag=new Float32Array(t);for(let r=0;r<t;r++){const s=-Math.PI*(r*r)/t;this._chirpReal[r]=Math.cos(s),this._chirpImag[r]=Math.sin(s)}const e=new Float32Array(this._m),i=new Float32Array(this._m);for(let r=0;r<t;r++)e[r]=this._chirpReal[r],i[r]=-this._chirpImag[r];for(let r=1;r<t;r++)e[this._m-r]=e[r],i[this._m-r]=i[r];this._internalFFT.transform(e,i),this._bReal=new Float32Array(this._internalFFT._real),this._bImag=new Float32Array(this._internalFFT._imag)}initHanning(){if(this._hanning)return;const t=this.size;this._hanning=new Float32Array(t);let e=0;for(let i=0;i<t;i++){const r=.5*(1-Math.cos(2*Math.PI*i/(t-1)));this._hanning[i]=r,e+=r*r}this._windowSumSq=e}transform(t,e){this.isPowerOfTwo?this.transformRadix2(t,e):this.transformBluestein(t,e)}transformRadix2(t,e){const i=this.size,r=this._rev,s=this._real,o=this._imag;if(t===s)for(let a=0;a<i;a++){const l=r[a];if(a<l){const u=s[a],h=o[a];s[a]=s[l],o[a]=o[l],s[l]=u,o[l]=h}}else for(let a=0;a<i;a++){const l=r[a];s[a]=t[l],o[a]=e?e[l]:0}for(let a=2;a<=i;a*=2){const l=a/2,u=-2*Math.PI/a,h=Math.cos(u),c=Math.sin(u);for(let g=0;g<i;g+=a){let d=1,f=0;for(let p=0;p<l;p++){const m=g+p,y=g+p+l,x=d*s[y]-f*o[y],b=d*o[y]+f*s[y],_=s[m],k=o[m];s[m]=_+x,o[m]=k+b,s[y]=_-x,o[y]=k-b;const v=d*h-f*c,w=d*c+f*h;d=v,f=w}}}}transformBluestein(t,e){const i=this.size,r=this._m,s=this._internalFFT,o=s._real,a=s._imag;o.fill(0),a.fill(0);for(let c=0;c<i;c++){const g=t[c],d=e?e[c]:0,f=this._chirpReal[c],p=this._chirpImag[c];o[c]=g*f-d*p,a[c]=g*p+d*f}s.transformRadix2(o,a);for(let c=0;c<r;c++){const g=s._real[c],d=s._imag[c],f=this._bReal[c],p=this._bImag[c];s._real[c]=g*f-d*p,s._imag[c]=g*p+d*f}const l=s._real,u=s._imag;for(let c=0;c<r;c++)u[c]=-u[c];s.transformRadix2(l,u);const h=1/r;for(let c=0;c<i;c++){const g=s._real[c]*h,d=-s._imag[c]*h,f=this._chirpReal[c],p=this._chirpImag[c];this._real[c]=g*f-d*p,this._imag[c]=g*p+d*f}}calculateSpectrum(t,e,i=!1){const r=this.size;let s=0;for(let h=0;h<r;h++)s+=t[h];const o=s/r;this._scratch||(this._scratch=new Float32Array(r));const a=this._scratch;if(i){this.initHanning();const h=this._hanning;for(let c=0;c<r;c++)a[c]=(t[c]-o)*h[c]}else for(let h=0;h<r;h++)a[h]=t[h]-o;this.transform(a);const l=e.length;let u=1/r;i&&this._windowSumSq>0&&(u=1/this._windowSumSq);for(let h=0;h<l;h++){const c=this._real[h],g=this._imag[h];e[h]+=(c*c+g*g)*u}}calculateSpectrumWindow(t,e,i,r=!1){const s=this.size;let o=0;for(let c=0;c<s;c++)o+=t[e+c];const a=o/s;this._scratch||(this._scratch=new Float32Array(s));const l=this._scratch;if(r){this.initHanning();const c=this._hanning;for(let g=0;g<s;g++)l[g]=(t[e+g]-a)*c[g]}else for(let c=0;c<s;c++)l[c]=t[e+c]-a;this.transform(l);const u=i.length;let h=1/s;r&&this._windowSumSq>0&&(h=1/this._windowSumSq);for(let c=0;c<u;c++){const g=this._real[c],d=this._imag[c];i[c]+=(g*g+d*d)*h}}}const qs={"Sony ILCE-7RM5":"0.82 -0.2976 -0.0719 -0.4296 1.2053 0.2532 -0.0429 0.1282 0.5774"};let Ai=null;async function Ks(n){return Ai||(Ai=(async()=>{if(typeof window.loadPyodide!="function")throw new Error("Pyodide missing: window.loadPyodide not found.");const t=await window.loadPyodide();return await t.loadPackage("numpy"),t})()),Ai}let Ti=null,xr=!1;async function Qi(){var e;Ti||(Ti=(async()=>{const i=await import("./joraw2-Bb3_vNP4.js");if(typeof i.default!="function")throw new Error("JoRaw2 WASM import failed");const r=new URL("/assets/joraw2-3YkywkGx.wasm",import.meta.url).href;return i.default({locateFile(s,o){return s.endsWith("joraw2.wasm")?r:o+s}})})());const t=(await Ti).LibRaw;if(!t)throw new Error("JoRaw2 class not found");if(!xr){const i=new t;try{if(typeof i.runtimeInfo!="function")throw new Error("JoRaw2 runtime identity is missing");const r=i.runtimeInfo();if((r==null?void 0:r.wrapper)!=="joraw2"||!String((r==null?void 0:r.librawVersion)||"").startsWith("0.22.2")||!(r!=null&&r.nikonHe)||!(r!=null&&r.nikonHeStar))throw new Error(`Unexpected JoRaw2 runtime: ${JSON.stringify(r)}`);xr=!0}finally{(e=i.delete)==null||e.call(i)}}return t}var $s=`#!/usr/bin/env python3
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
`,Js=`#!/usr/bin/env python3
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
`,Zs=`#!/usr/bin/env python3
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
`;const br=512,ta=qs["Sony ILCE-7RM5"].split(/\s+/).map(Number).filter(Number.isFinite);let Ii=null;function ea(n){if(n.byteLength<8)return null;const t=n.getUint16(0,!1);return t===18761?!0:t===19789?!1:null}function Zr(n,t,e){const i=Math.min(n.length,t+e);let r="";for(let s=t;s<i;s++){const o=n[s];if(o===0)break;r+=String.fromCharCode(o)}return r.trim()}function Qt(n,t,e,i,r){const s=e===1||e===2||e===7?1:e===3||e===8?2:e===4||e===9?4:0;if(!s)return[];const o=s*i,a=o<=4?r:n.getUint32(r,t);if(a<0||a+o>n.byteLength)return[];const l=[];for(let u=0;u<i;u++){const h=a+u*s;e===1||e===2||e===7?l.push(n.getUint8(h)):e===3?l.push(n.getUint16(h,t)):e===8?l.push(n.getInt16(h,t)):e===4?l.push(n.getUint32(h,t)):e===9&&l.push(n.getInt32(h,t))}return l}function _r(n,t,e,i,r,s){if(i!==2||r<=0)return"";const o=r<=4?s:t.getUint32(s,e);return o<0||o>=n.length?"":Zr(n,o,r)}function ts(n){const t=new Uint8Array(n),e=new DataView(n),i=ea(e);if(i===null||e.getUint16(2,i)!==42)return null;const s=c=>e.getUint16(c,i),o=c=>e.getUint32(c,i),a=[o(4)],l=new Set;let u="",h="";for(;a.length;){const c=a.pop();if(l.has(c)||c<=0||c+2>e.byteLength)continue;l.add(c);const g=s(c);if(c+2+g*12+4>e.byteLength)continue;const d=new Map;for(let _=0;_<g;_++){const k=c+2+_*12,v=s(k),w=s(k+2),P=o(k+4),A=k+8;d.set(v,{type:w,count:P,valueOffset:A})}const f=d.get(271),p=d.get(272);f&&!u&&(u=_r(t,e,i,f.type,f.count,f.valueOffset)),p&&!h&&(h=_r(t,e,i,p.type,p.count,p.valueOffset));const m=d.get(330);if(m){const _=Qt(e,i,m.type,m.count,m.valueOffset);for(const k of _)a.push(k)}const y=d.get(259),x=d.get(262);if(y&&x){const _=Qt(e,i,y.type,y.count,y.valueOffset)[0],k=Qt(e,i,x.type,x.count,x.valueOffset)[0];if(_===32766&&k===32803){const v=Qt(e,i,d.get(256).type,d.get(256).count,d.get(256).valueOffset)[0],w=Qt(e,i,d.get(257).type,d.get(257).count,d.get(257).valueOffset)[0],P=Qt(e,i,d.get(258).type,d.get(258).count,d.get(258).valueOffset)[0],A=Qt(e,i,d.get(273).type,d.get(273).count,d.get(273).valueOffset)[0],F=Qt(e,i,d.get(279).type,d.get(279).count,d.get(279).valueOffset)[0],C=d.get(33422)?Qt(e,i,d.get(33422).type,d.get(33422).count,d.get(33422).valueOffset):[0,1,1,2];d.get(29456)&&Qt(e,i,d.get(29456).type,d.get(29456).count,d.get(29456).valueOffset);const M=d.get(50717)?Qt(e,i,d.get(50717).type,d.get(50717).count,d.get(50717).valueOffset)[0]:16383,S=d.get(50719)?Qt(e,i,d.get(50719).type,d.get(50719).count,d.get(50719).valueOffset):[],I=d.get(50720)?Qt(e,i,d.get(50720).type,d.get(50720).count,d.get(50720).valueOffset):[];if(A+br+16>t.length||A+F>t.length)return null;const N=A+br,E=Zr(t,N,4),X=t[N+8]<<8|t[N+9],U=t[N+10]<<8|t[N+11],L=t[N+12]<<8|t[N+13],B=t[N+14]<<8|t[N+15],O=L>>4&63,z=B>>13,J=B>>10&3,V=U*2,q=F>=4?(t[A]|t[A+1]<<8|t[A+2]<<16|t[A+3]<<24)>>>0:0,W=X===v&&V===w;let tt=!1;if(q>=1&&q<=16&&F>=8+q*24){const H=new Map,T=new Map;let Q=!0;for(let G=0;G<q;G++){const Y=A+8+G*24,et=e.getUint32(Y+8,!0),nt=e.getUint32(Y+12,!0),it=e.getUint32(Y+16,!0),ct=e.getUint32(Y+20,!0);if(!it||!ct||et+it>v||nt+ct>w){Q=!1;break}const dt=H.get(nt);if(dt!==void 0&&dt!==ct){Q=!1;break}H.set(nt,ct),T.set(nt,(T.get(nt)||0)+it)}if(Q){const G=Array.from(H.keys()).sort((et,nt)=>et-nt);let Y=0;for(const et of G){if(et!==Y||T.get(et)!==v){Q=!1;break}Y+=H.get(et)}tt=Q&&Y===w}}const at=q>=1&&q<=16&&X>0&&V>0&&v%X===0&&w%V===0&&q===v/X*(w/V);if(E!=="A000"&&E!=="0000"||!W&&!at&&!tt||O!==16||z!==3||J!==3)return null;const D=[1024,1024,1024,1024],Z=C.length>=4?C.slice(0,4).map(H=>H===0?"R":H===2?"B":H===1?"G":"?").join(""):"";return{width:v,height:w,bitsPerSample:P,compression:_,photometric:k,blackLevel:D,whiteLevel:Number(M||16383),cfaPattern:Z,defaultCropOrigin:S.length>=2?[Number(S[0]),Number(S[1])]:void 0,defaultCropSize:I.length>=2?[Number(I[0]),Number(I[1])]:void 0,make:u||"SONY",model:h||"ILCE-7M5"}}}const b=o(c+2+g*12);b&&a.push(b)}return null}async function na(n){return Ii||(Ii=(async()=>{const t=await Ks();return t.__jtrSonyCrawHqDecoderReady||(await t.FS.mkdirTree("/sony_craw_hq"),await t.FS.writeFile("/sony_craw_hq/llvc3_bitstream_probe.py",$s),await t.FS.writeFile("/sony_craw_hq/llvc3_entropy.py",Js),await t.FS.writeFile("/sony_craw_hq/llvc3_math.py",Zs),await t.runPythonAsync(`
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
`),t.__jtrSonyCrawHqDecoderReady=!0),t})()),Ii}function ia(n){return ts(n)}function ra(n){return n==="ILCE-7M5"?ta:null}function es(n,t,e,i){var l;if(n.length!==t.width*t.height)throw new Error(`Sony cRAW HQ decoded size mismatch: got ${n.length}, expected ${t.width*t.height}`);const r=t.model||"ILCE-7M5",s=r.startsWith("Sony ")?r:`Sony ${r}`,o=ra(r),a={...i||{},make:t.make||(i==null?void 0:i.camera_make)||"SONY",model:r,camera_make:t.make||(i==null?void 0:i.camera_make)||"SONY",camera_model:r,UniqueCameraModel:s,sourceFormat:e==="joraw2-wasm"?"Sony cRAW HQ / LLVC3 (JoRaw2 WASM)":"Sony cRAW HQ / LLVC3",sonyCrawHq:{...t,decodeBackend:e},color_desc:t.cfaPattern,black_level_per_channel:t.blackLevel,white_level:t.whiteLevel,color_matrix:o&&o.length===9?o:void 0,idata:{filters:2492765332,colors:3},color_data:{...(i==null?void 0:i.color_data)||{},black:1024,cblack_rawpy_style:t.blackLevel,dng_levels:{...((l=i==null?void 0:i.color_data)==null?void 0:l.dng_levels)||{},dng_cblack:t.blackLevel,dng_whitelevel:t.whiteLevel}}};return{data:n,width:t.width,height:t.height,bayerPattern:t.cfaPattern,blackLevels:t.blackLevel,whiteLevel:t.whiteLevel,metadata:a,isThreePlane:!1,isXTrans:!1}}async function sa(n,t,e){const i=typeof performance<"u"?performance.now():Date.now(),r=await Qi(),s=typeof performance<"u"?performance.now():Date.now(),o=new r;try{const a=new Uint8Array(n);await o.open(a,{});const l=typeof performance<"u"?performance.now():Date.now();let u=null;try{u=await o.metadata(!0)}catch(m){console.warn("[Sony cRAW HQ] fast WASM metadata read failed",m)}const h=typeof performance<"u"?performance.now():Date.now(),c=o.getRawImage(),g=typeof performance<"u"?performance.now():Date.now();if(!c||!c.data)throw new Error("Sony cRAW HQ LibRaw WASM returned no raw image");const d=c.data instanceof Uint16Array?c.data:new Uint16Array(c.data.buffer,c.data.byteOffset||0,c.data.byteLength/2),f=typeof performance<"u"?performance.now():Date.now();if(c.width!==t.width||c.height!==t.height)throw new Error(`Sony cRAW HQ LibRaw WASM dimensions mismatch: got ${c.width}x${c.height}, expected ${t.width}x${t.height}`);const p=typeof performance<"u"?performance.now():Date.now();return console.info("[Sony cRAW HQ] fast decode timings",{width:t.width,height:t.height,backend:"joraw2-wasm",wasmReadyMs:Math.round(s-i),openMs:Math.round(l-s),metadataMs:Math.round(h-l),unpackMs:Math.round(g-h),copyMs:Math.round(f-g),totalMs:Math.round(p-i)}),{rawImageData:es(d,t,"joraw2-wasm",u),info:t}}finally{typeof o.delete=="function"?o.delete():typeof o.close=="function"&&o.close()}}async function aa(n,t,e){const i=typeof performance<"u"?performance.now():Date.now(),r=await na(),s=typeof performance<"u"?performance.now():Date.now(),o=new Uint8Array(n),a=await fetch(new URL("/assets/sony_llvc3_static_lut4096_padded_u16-FsVBk-IV.bin",import.meta.url));if(!a.ok)throw new Error(`Failed to load Sony LLVC3 sample LUT: HTTP ${a.status}`);const l=new Uint8Array(await a.arrayBuffer()),u=typeof performance<"u"?performance.now():Date.now();r.globals.set("jtr_sony_arw_bytes",o),r.globals.set("jtr_sony_lut_bytes",l);const h=await r.runPythonAsync("jtr_decode_sony_craw_hq(jtr_sony_arw_bytes.to_py(), jtr_sony_lut_bytes.to_py())"),c=typeof performance<"u"?performance.now():Date.now(),g=h.toJs();typeof h.destroy=="function"&&h.destroy(),r.globals.delete("jtr_sony_arw_bytes"),r.globals.delete("jtr_sony_lut_bytes");const d=new Uint8Array(g.byteLength);d.set(g);const f=new Uint16Array(d.buffer),p=typeof performance<"u"?performance.now():Date.now(),m=typeof performance<"u"?performance.now():Date.now();return console.info("[Sony cRAW HQ] decode timings",{width:t.width,height:t.height,backend:"pyodide",pyodideReadyMs:Math.round(s-i),lutLoadMs:Math.round(u-s),llvc3DecodeMs:Math.round(c-u),copyMs:Math.round(p-c),totalMs:Math.round(m-i)}),{rawImageData:es(f,t,"pyodide"),info:t}}async function oa(n,t){const e=ts(n);if(!e)return null;try{return await sa(n,e,t)}catch(i){return console.warn("[Sony cRAW HQ] fast WASM decode failed; falling back to Pyodide",i),aa(n,e)}}async function la(n,t){return oa(n,t)}const ca=["RGGB","BGGR","GRBG","GBRG"],ua=new Set(ca);function ha(n){if(!n||typeof n!="object"||Array.isArray(n)||ArrayBuffer.isView(n))return n;const t=n;return t.value??t.values??t.description??n}function Je(n){const t=ha(n);let e="";if(typeof t=="string")e=t.toUpperCase().replace(/[^RGB012]/g,"");else if(typeof t=="number")e=String(t);else if(Array.isArray(t)||ArrayBuffer.isView(t))e=Array.from(t).map(String).join("");else return null;return/^[012]{4}$/.test(e)&&(e=e.replace(/0/g,"R").replace(/1/g,"G").replace(/2/g,"B")),ua.has(e)?e:null}function Mr(n){var t;return n?[n.cfa_pattern,n.cfaPattern,n.BayerPattern,n.CFAPattern2,n.CFAPattern,(t=n.idata)==null?void 0:t.cfa_pattern]:[]}function da(n,t,e){const i=Je(n.bayerPattern),r=n.bayerPatternSource==="manual";if(i&&!r)return{pattern:i,source:n.bayerPatternSource||"decoder"};if(!r)for(const o of Mr(n.metadata)){const a=Je(o);if(a)return{pattern:a,source:"libraw-metadata"}}for(const o of Mr(t)){const a=Je(o);if(a)return{pattern:a,source:"metadata"}}const s=Je(e)||(r?i:null);return s?{pattern:s,source:"manual"}:{pattern:null,source:null}}function fa(n,t){return n.bayerPattern=t.pattern||"",n.bayerPatternSource=t.source||void 0,t.pattern&&t.source!=="manual"&&(n.metadata={...n.metadata||{},cfa_pattern:t.pattern,cfaPattern:t.pattern}),t}function ns(n,t){return(t&1)<<1|n&1}function pa(n,t,e){const i=Je(n);if(!i)throw new Error("Bayer CFA pattern is unresolved.");return i[ns(t,e)]}const xe=n=>{const t=Number(n);return Number.isFinite(t)?Math.max(0,t):0};function $e(n){if(!n||typeof n.length!="number")return null;const t=Array.from(n);return t.length<4?null:[xe(t[0]),xe(t[1]),xe(t[2]),xe(t[3])]}function wr(n){if(!n||n.source!=="libraw")return null;const t=$e(n.channelOffsets),e=$e(n.channelLevels),i=$e(n.siteColorIndices),r=$e(n.siteBaseLevels),s=$e(n.siteLevels);if(!t||!e||!i||!r||!s)return null;const o=Math.max(0,Math.floor(xe(n.repeatRows))),a=Math.max(0,Math.floor(xe(n.repeatCols))),l=o*a,u=n.repeatValues&&typeof n.repeatValues.length=="number"?Array.from(n.repeatValues,xe):[],h=l>0&&l<=4098&&u.length>=l;return{source:"libraw",common:xe(n.common),channelOffsets:t,channelLevels:e,siteColorIndices:i,siteBaseLevels:r,siteLevels:s,repeatRows:h?o:0,repeatCols:h?a:0,repeatOriginY:Math.max(0,Math.floor(xe(n.repeatOriginY))),repeatOriginX:Math.max(0,Math.floor(xe(n.repeatOriginX))),repeatValues:h?u.slice(0,l):[]}}var oi=typeof self<"u"?self:global;const Tn=typeof navigator<"u",ma=Tn&&typeof HTMLImageElement>"u",Zn=!(typeof global>"u"||typeof process>"u"||!process.versions||!process.versions.node),li=oi.Buffer,Yn=oi.BigInt,ci=!!li,ga=n=>n;function ti(n,t=ga){if(Zn)try{return typeof require=="function"?Promise.resolve(t(require(n))):import(n).then(t)}catch{console.warn(`Couldn't load ${n}`)}}let qi=oi.fetch;const ya=n=>qi=n;if(!oi.fetch){const n=ti("http",i=>i),t=ti("https",i=>i),e=(i,{headers:r}={})=>new Promise(async(s,o)=>{let{port:a,hostname:l,pathname:u,protocol:h,search:c}=new URL(i);const g={method:"GET",hostname:l,path:encodeURI(u)+c,headers:r};a!==""&&(g.port=Number(a));const d=(h==="https:"?await t:await n).request(g,f=>{if(f.statusCode===301||f.statusCode===302){let p=new URL(f.headers.location,i).toString();return e(p,{headers:r}).then(s).catch(o)}s({status:f.statusCode,arrayBuffer:()=>new Promise(p=>{let m=[];f.on("data",y=>m.push(y)),f.on("end",()=>p(Buffer.concat(m)))})})});d.on("error",o),d.end()});ya(e)}function st(n,t,e){return t in n?Object.defineProperty(n,t,{value:e,enumerable:!0,configurable:!0,writable:!0}):n[t]=e,n}const ei=n=>is(n)?void 0:n,xa=n=>n!==void 0;function is(n){return n===void 0||(n instanceof Map?n.size===0:Object.values(n).filter(xa).length===0)}function kt(n){let t=new Error(n);throw delete t.stack,t}function Ze(n){return(n=function(t){for(;t.endsWith("\0");)t=t.slice(0,-1);return t}(n).trim())===""?void 0:n}function Vi(n){let t=function(e){let i=0;return e.ifd0.enabled&&(i+=1024),e.exif.enabled&&(i+=2048),e.makerNote&&(i+=2048),e.userComment&&(i+=1024),e.gps.enabled&&(i+=512),e.interop.enabled&&(i+=100),e.ifd1.enabled&&(i+=1024),i+2048}(n);return n.jfif.enabled&&(t+=50),n.xmp.enabled&&(t+=2e4),n.iptc.enabled&&(t+=14e3),n.icc.enabled&&(t+=6e3),t}const Xi=n=>String.fromCharCode.apply(null,n),Sr=typeof TextDecoder<"u"?new TextDecoder("utf-8"):void 0;function rs(n){return Sr?Sr.decode(n):ci?Buffer.from(n).toString("utf8"):decodeURIComponent(escape(Xi(n)))}class jt{static from(t,e){return t instanceof this&&t.le===e?t:new jt(t,void 0,void 0,e)}constructor(t,e=0,i,r){if(typeof r=="boolean"&&(this.le=r),Array.isArray(t)&&(t=new Uint8Array(t)),t===0)this.byteOffset=0,this.byteLength=0;else if(t instanceof ArrayBuffer){i===void 0&&(i=t.byteLength-e);let s=new DataView(t,e,i);this._swapDataView(s)}else if(t instanceof Uint8Array||t instanceof DataView||t instanceof jt){i===void 0&&(i=t.byteLength-e),(e+=t.byteOffset)+i>t.byteOffset+t.byteLength&&kt("Creating view outside of available memory in ArrayBuffer");let s=new DataView(t.buffer,e,i);this._swapDataView(s)}else if(typeof t=="number"){let s=new DataView(new ArrayBuffer(t));this._swapDataView(s)}else kt("Invalid input argument for BufferView: "+t)}_swapArrayBuffer(t){this._swapDataView(new DataView(t))}_swapBuffer(t){this._swapDataView(new DataView(t.buffer,t.byteOffset,t.byteLength))}_swapDataView(t){this.dataView=t,this.buffer=t.buffer,this.byteOffset=t.byteOffset,this.byteLength=t.byteLength}_lengthToEnd(t){return this.byteLength-t}set(t,e,i=jt){return t instanceof DataView||t instanceof jt?t=new Uint8Array(t.buffer,t.byteOffset,t.byteLength):t instanceof ArrayBuffer&&(t=new Uint8Array(t)),t instanceof Uint8Array||kt("BufferView.set(): Invalid data argument."),this.toUint8().set(t,e),new i(this,e,t.byteLength)}subarray(t,e){return e=e||this._lengthToEnd(t),new jt(this,t,e)}toUint8(){return new Uint8Array(this.buffer,this.byteOffset,this.byteLength)}getUint8Array(t,e){return new Uint8Array(this.buffer,this.byteOffset+t,e)}getString(t=0,e=this.byteLength){return rs(this.getUint8Array(t,e))}getLatin1String(t=0,e=this.byteLength){let i=this.getUint8Array(t,e);return Xi(i)}getUnicodeString(t=0,e=this.byteLength){const i=[];for(let r=0;r<e&&t+r<this.byteLength;r+=2)i.push(this.getUint16(t+r));return Xi(i)}getInt8(t){return this.dataView.getInt8(t)}getUint8(t){return this.dataView.getUint8(t)}getInt16(t,e=this.le){return this.dataView.getInt16(t,e)}getInt32(t,e=this.le){return this.dataView.getInt32(t,e)}getUint16(t,e=this.le){return this.dataView.getUint16(t,e)}getUint32(t,e=this.le){return this.dataView.getUint32(t,e)}getFloat32(t,e=this.le){return this.dataView.getFloat32(t,e)}getFloat64(t,e=this.le){return this.dataView.getFloat64(t,e)}getFloat(t,e=this.le){return this.dataView.getFloat32(t,e)}getDouble(t,e=this.le){return this.dataView.getFloat64(t,e)}getUintBytes(t,e,i){switch(e){case 1:return this.getUint8(t,i);case 2:return this.getUint16(t,i);case 4:return this.getUint32(t,i);case 8:return this.getUint64&&this.getUint64(t,i)}}getUint(t,e,i){switch(e){case 8:return this.getUint8(t,i);case 16:return this.getUint16(t,i);case 32:return this.getUint32(t,i);case 64:return this.getUint64&&this.getUint64(t,i)}}toString(t){return this.dataView.toString(t,this.constructor.name)}ensureChunk(){}}function Gi(n,t){kt(`${n} '${t}' was not loaded, try using full build of exifr.`)}class Ki extends Map{constructor(t){super(),this.kind=t}get(t,e){return this.has(t)||Gi(this.kind,t),e&&(t in e||function(i,r){kt(`Unknown ${i} '${r}'.`)}(this.kind,t),e[t].enabled||Gi(this.kind,t)),super.get(t)}keyList(){return Array.from(this.keys())}}var be=new Ki("file parser"),Ct=new Ki("segment parser"),Me=new Ki("file reader");function ba(n,t){return typeof n=="string"?vr(n,t):Tn&&!ma&&n instanceof HTMLImageElement?vr(n.src,t):n instanceof Uint8Array||n instanceof ArrayBuffer||n instanceof DataView?new jt(n):Tn&&n instanceof Blob?zi(n,t,"blob",ln):void kt("Invalid input argument")}function vr(n,t){return(e=n).startsWith("data:")||e.length>1e4?Yi(n,t,"base64"):Zn&&n.includes("://")?zi(n,t,"url",on):Zn?Yi(n,t,"fs"):Tn?zi(n,t,"url",on):void kt("Invalid input argument");var e}async function zi(n,t,e,i){return Me.has(e)?Yi(n,t,e):i?async function(r,s){let o=await s(r);return new jt(o)}(n,i):void kt(`Parser ${e} is not loaded`)}async function Yi(n,t,e){let i=new(Me.get(e))(n,t);return await i.read(),i}const on=n=>qi(n).then(t=>t.arrayBuffer()),ln=n=>new Promise((t,e)=>{let i=new FileReader;i.onloadend=()=>t(i.result||new ArrayBuffer),i.onerror=e,i.readAsArrayBuffer(n)});class _a extends Map{get tagKeys(){return this.allKeys||(this.allKeys=Array.from(this.keys())),this.allKeys}get tagValues(){return this.allValues||(this.allValues=Array.from(this.values())),this.allValues}}function Pt(n,t,e){let i=new _a;for(let[r,s]of e)i.set(r,s);if(Array.isArray(t))for(let r of t)n.set(r,i);else n.set(t,i);return i}function cn(n,t,e){let i,r=n.get(t);for(i of e)r.set(i[0],i[1])}const Rt=new Map,le=new Map,Be=new Map,Re=["chunked","firstChunkSize","firstChunkSizeNode","firstChunkSizeBrowser","chunkSize","chunkLimit"],bn=["jfif","xmp","icc","iptc","ihdr"],un=["tiff",...bn],St=["ifd0","ifd1","exif","gps","interop"],Ne=[...un,...St],Le=["makerNote","userComment"],_n=["translateKeys","translateValues","reviveValues","multiSegment"],Ee=[..._n,"sanitize","mergeOutput","silentErrors"];class ss{get translate(){return this.translateKeys||this.translateValues||this.reviveValues}}class Cn extends ss{get needed(){return this.enabled||this.deps.size>0}constructor(t,e,i,r){if(super(),st(this,"enabled",!1),st(this,"skip",new Set),st(this,"pick",new Set),st(this,"deps",new Set),st(this,"translateKeys",!1),st(this,"translateValues",!1),st(this,"reviveValues",!1),this.key=t,this.enabled=e,this.parse=this.enabled,this.applyInheritables(r),this.canBeFiltered=St.includes(t),this.canBeFiltered&&(this.dict=Rt.get(t)),i!==void 0)if(Array.isArray(i))this.parse=this.enabled=!0,this.canBeFiltered&&i.length>0&&this.translateTagSet(i,this.pick);else if(typeof i=="object"){if(this.enabled=!0,this.parse=i.parse!==!1,this.canBeFiltered){let{pick:s,skip:o}=i;s&&s.length>0&&this.translateTagSet(s,this.pick),o&&o.length>0&&this.translateTagSet(o,this.skip)}this.applyInheritables(i)}else i===!0||i===!1?this.parse=this.enabled=i:kt(`Invalid options argument: ${i}`)}applyInheritables(t){let e,i;for(e of _n)i=t[e],i!==void 0&&(this[e]=i)}translateTagSet(t,e){if(this.dict){let i,r,{tagKeys:s,tagValues:o}=this.dict;for(i of t)typeof i=="string"?(r=o.indexOf(i),r===-1&&(r=s.indexOf(Number(i))),r!==-1&&e.add(Number(s[r]))):e.add(i)}else for(let i of t)e.add(i)}finalizeFilters(){!this.enabled&&this.deps.size>0?(this.enabled=!0,ni(this.pick,this.deps)):this.enabled&&this.pick.size>0&&ni(this.pick,this.deps)}}var Ot={jfif:!1,tiff:!0,xmp:!1,icc:!1,iptc:!1,ifd0:!0,ifd1:!1,exif:!0,gps:!0,interop:!1,ihdr:void 0,makerNote:!1,userComment:!1,multiSegment:!1,skip:[],pick:[],translateKeys:!0,translateValues:!0,reviveValues:!0,sanitize:!0,mergeOutput:!0,silentErrors:!0,chunked:!0,firstChunkSize:void 0,firstChunkSizeNode:512,firstChunkSizeBrowser:65536,chunkSize:65536,chunkLimit:5},Pr=new Map;class hn extends ss{static useCached(t){let e=Pr.get(t);return e!==void 0||(e=new this(t),Pr.set(t,e)),e}constructor(t){super(),t===!0?this.setupFromTrue():t===void 0?this.setupFromUndefined():Array.isArray(t)?this.setupFromArray(t):typeof t=="object"?this.setupFromObject(t):kt(`Invalid options argument ${t}`),this.firstChunkSize===void 0&&(this.firstChunkSize=Tn?this.firstChunkSizeBrowser:this.firstChunkSizeNode),this.mergeOutput&&(this.ifd1.enabled=!1),this.filterNestedSegmentTags(),this.traverseTiffDependencyTree(),this.checkLoadedPlugins()}setupFromUndefined(){let t;for(t of Re)this[t]=Ot[t];for(t of Ee)this[t]=Ot[t];for(t of Le)this[t]=Ot[t];for(t of Ne)this[t]=new Cn(t,Ot[t],void 0,this)}setupFromTrue(){let t;for(t of Re)this[t]=Ot[t];for(t of Ee)this[t]=Ot[t];for(t of Le)this[t]=!0;for(t of Ne)this[t]=new Cn(t,!0,void 0,this)}setupFromArray(t){let e;for(e of Re)this[e]=Ot[e];for(e of Ee)this[e]=Ot[e];for(e of Le)this[e]=Ot[e];for(e of Ne)this[e]=new Cn(e,!1,void 0,this);this.setupGlobalFilters(t,void 0,St)}setupFromObject(t){let e;for(e of(St.ifd0=St.ifd0||St.image,St.ifd1=St.ifd1||St.thumbnail,Object.assign(this,t),Re))this[e]=Ri(t[e],Ot[e]);for(e of Ee)this[e]=Ri(t[e],Ot[e]);for(e of Le)this[e]=Ri(t[e],Ot[e]);for(e of un)this[e]=new Cn(e,Ot[e],t[e],this);for(e of St)this[e]=new Cn(e,Ot[e],t[e],this.tiff);this.setupGlobalFilters(t.pick,t.skip,St,Ne),t.tiff===!0?this.batchEnableWithBool(St,!0):t.tiff===!1?this.batchEnableWithUserValue(St,t):Array.isArray(t.tiff)?this.setupGlobalFilters(t.tiff,void 0,St):typeof t.tiff=="object"&&this.setupGlobalFilters(t.tiff.pick,t.tiff.skip,St)}batchEnableWithBool(t,e){for(let i of t)this[i].enabled=e}batchEnableWithUserValue(t,e){for(let i of t){let r=e[i];this[i].enabled=r!==!1&&r!==void 0}}setupGlobalFilters(t,e,i,r=i){if(t&&t.length){for(let o of r)this[o].enabled=!1;let s=Cr(t,i);for(let[o,a]of s)ni(this[o].pick,a),this[o].enabled=!0}else if(e&&e.length){let s=Cr(e,i);for(let[o,a]of s)ni(this[o].skip,a)}}filterNestedSegmentTags(){let{ifd0:t,exif:e,xmp:i,iptc:r,icc:s}=this;this.makerNote?e.deps.add(37500):e.skip.add(37500),this.userComment?e.deps.add(37510):e.skip.add(37510),i.enabled||t.skip.add(700),r.enabled||t.skip.add(33723),s.enabled||t.skip.add(34675)}traverseTiffDependencyTree(){let{ifd0:t,exif:e,gps:i,interop:r}=this;r.needed&&(e.deps.add(40965),t.deps.add(40965)),e.needed&&t.deps.add(34665),i.needed&&t.deps.add(34853),this.tiff.enabled=St.some(s=>this[s].enabled===!0)||this.makerNote||this.userComment;for(let s of St)this[s].finalizeFilters()}get onlyTiff(){return!bn.map(t=>this[t].enabled).some(t=>t===!0)&&this.tiff.enabled}checkLoadedPlugins(){for(let t of un)this[t].enabled&&!Ct.has(t)&&Gi("segment parser",t)}}function Cr(n,t){let e,i,r,s,o=[];for(r of t){for(s of(e=Rt.get(r),i=[],e))(n.includes(s[0])||n.includes(s[1]))&&i.push(s[0]);i.length&&o.push([r,i])}return o}function Ri(n,t){return n!==void 0?n:t!==void 0?t:void 0}function ni(n,t){for(let e of t)n.add(e)}st(hn,"default",Ot);class De{constructor(t){st(this,"parsers",{}),st(this,"output",{}),st(this,"errors",[]),st(this,"pushToErrors",e=>this.errors.push(e)),this.options=hn.useCached(t)}async read(t){this.file=await ba(t,this.options)}setup(){if(this.fileParser)return;let{file:t}=this,e=t.getUint16(0);for(let[i,r]of be)if(r.canHandle(t,e))return this.fileParser=new r(this.options,this.file,this.parsers),t[i]=!0;this.file.close&&this.file.close(),kt("Unknown file format")}async parse(){let{output:t,errors:e}=this;return this.setup(),this.options.silentErrors?(await this.executeParsers().catch(this.pushToErrors),e.push(...this.fileParser.errors)):await this.executeParsers(),this.file.close&&this.file.close(),this.options.silentErrors&&e.length>0&&(t.errors=e),ei(t)}async executeParsers(){let{output:t}=this;await this.fileParser.parse();let e=Object.values(this.parsers).map(async i=>{let r=await i.parse();i.assignToOutput(t,r)});this.options.silentErrors&&(e=e.map(i=>i.catch(this.pushToErrors))),await Promise.all(e)}async extractThumbnail(){this.setup();let{options:t,file:e}=this,i=Ct.get("tiff",t);var r;if(e.tiff?r={start:0,type:"tiff"}:e.jpeg&&(r=await this.fileParser.getOrFindSegment("tiff")),r===void 0)return;let s=await this.fileParser.ensureSegmentChunk(r),o=this.parsers.tiff=new i(s,t,e),a=await o.extractThumbnail();return e.close&&e.close(),a}}async function ui(n,t){let e=new De(t);return await e.read(n),e.parse()}var Ma=Object.freeze({__proto__:null,parse:ui,Exifr:De,fileParsers:be,segmentParsers:Ct,fileReaders:Me,tagKeys:Rt,tagValues:le,tagRevivers:Be,createDictionary:Pt,extendDictionary:cn,fetchUrlAsArrayBuffer:on,readBlobAsArrayBuffer:ln,chunkedProps:Re,otherSegments:bn,segments:un,tiffBlocks:St,segmentsAndBlocks:Ne,tiffExtractables:Le,inheritables:_n,allFormatters:Ee,Options:hn});class hi{constructor(t,e,i){st(this,"errors",[]),st(this,"ensureSegmentChunk",async r=>{let s=r.start,o=r.size||65536;if(this.file.chunked)if(this.file.available(s,o))r.chunk=this.file.subarray(s,o);else try{r.chunk=await this.file.readChunk(s,o)}catch(a){kt(`Couldn't read segment: ${JSON.stringify(r)}. ${a.message}`)}else this.file.byteLength>s+o?r.chunk=this.file.subarray(s,o):r.size===void 0?r.chunk=this.file.subarray(s):kt("Segment unreachable: "+JSON.stringify(r));return r.chunk}),this.extendOptions&&this.extendOptions(t),this.options=t,this.file=e,this.parsers=i}injectSegment(t,e){this.options[t].enabled&&this.createParser(t,e)}createParser(t,e){let i=new(Ct.get(t))(e,this.options,this.file);return this.parsers[t]=i}createParsers(t){for(let e of t){let{type:i,chunk:r}=e,s=this.options[i];if(s&&s.enabled){let o=this.parsers[i];o&&o.append||o||this.createParser(i,r)}}}async readSegments(t){let e=t.map(this.ensureSegmentChunk);await Promise.all(e)}}class oe{static findPosition(t,e){let i=t.getUint16(e+2)+2,r=typeof this.headerLength=="function"?this.headerLength(t,e,i):this.headerLength,s=e+r,o=i-r;return{offset:e,length:i,headerLength:r,start:s,size:o,end:s+o}}static parse(t,e={}){return new this(t,new hn({[this.type]:e}),t).parse()}normalizeInput(t){return t instanceof jt?t:new jt(t)}constructor(t,e={},i){st(this,"errors",[]),st(this,"raw",new Map),st(this,"handleError",r=>{if(!this.options.silentErrors)throw r;this.errors.push(r.message)}),this.chunk=this.normalizeInput(t),this.file=i,this.type=this.constructor.type,this.globalOptions=this.options=e,this.localOptions=e[this.type],this.canTranslate=this.localOptions&&this.localOptions.translate}translate(){this.canTranslate&&(this.translated=this.translateBlock(this.raw,this.type))}get output(){return this.translated?this.translated:this.raw?Object.fromEntries(this.raw):void 0}translateBlock(t,e){let i=Be.get(e),r=le.get(e),s=Rt.get(e),o=this.options[e],a=o.reviveValues&&!!i,l=o.translateValues&&!!r,u=o.translateKeys&&!!s,h={};for(let[c,g]of t)a&&i.has(c)?g=i.get(c)(g):l&&r.has(c)&&(g=this.translateValue(g,r.get(c))),u&&s.has(c)&&(c=s.get(c)||c),h[c]=g;return h}translateValue(t,e){return e[t]||e.DEFAULT||t}assignToOutput(t,e){this.assignObjectToOutput(t,this.constructor.type,e)}assignObjectToOutput(t,e,i){if(this.globalOptions.mergeOutput)return Object.assign(t,i);t[e]?Object.assign(t[e],i):t[e]=i}}st(oe,"headerLength",4),st(oe,"type",void 0),st(oe,"multiSegment",!1),st(oe,"canHandle",()=>!1);function wa(n){return n===192||n===194||n===196||n===219||n===221||n===218||n===254}function Sa(n){return n>=224&&n<=239}function va(n,t,e){for(let[i,r]of Ct)if(r.canHandle(n,t,e))return i}class kr extends hi{constructor(...t){super(...t),st(this,"appSegments",[]),st(this,"jpegSegments",[]),st(this,"unknownSegments",[])}static canHandle(t,e){return e===65496}async parse(){await this.findAppSegments(),await this.readSegments(this.appSegments),this.mergeMultiSegments(),this.createParsers(this.mergedAppSegments||this.appSegments)}setupSegmentFinderArgs(t){t===!0?(this.findAll=!0,this.wanted=new Set(Ct.keyList())):(t=t===void 0?Ct.keyList().filter(e=>this.options[e].enabled):t.filter(e=>this.options[e].enabled&&Ct.has(e)),this.findAll=!1,this.remaining=new Set(t),this.wanted=new Set(t)),this.unfinishedMultiSegment=!1}async findAppSegments(t=0,e){this.setupSegmentFinderArgs(e);let{file:i,findAll:r,wanted:s,remaining:o}=this;if(!r&&this.file.chunked&&(r=Array.from(s).some(a=>{let l=Ct.get(a),u=this.options[a];return l.multiSegment&&u.multiSegment}),r&&await this.file.readWhole()),t=this.findAppSegmentsInRange(t,i.byteLength),!this.options.onlyTiff&&i.chunked){let a=!1;for(;o.size>0&&!a&&(i.canReadNextChunk||this.unfinishedMultiSegment);){let{nextChunkOffset:l}=i,u=this.appSegments.some(h=>!this.file.available(h.offset||h.start,h.length||h.size));if(a=t>l&&!u?!await i.readNextChunk(t):!await i.readNextChunk(l),(t=this.findAppSegmentsInRange(t,i.byteLength))===void 0)return}}}findAppSegmentsInRange(t,e){e-=2;let i,r,s,o,a,l,{file:u,findAll:h,wanted:c,remaining:g,options:d}=this;for(;t<e;t++)if(u.getUint8(t)===255){if(i=u.getUint8(t+1),Sa(i)){if(r=u.getUint16(t+2),s=va(u,t,r),s&&c.has(s)&&(o=Ct.get(s),a=o.findPosition(u,t),l=d[s],a.type=s,this.appSegments.push(a),!h&&(o.multiSegment&&l.multiSegment?(this.unfinishedMultiSegment=a.chunkNumber<a.chunkCount,this.unfinishedMultiSegment||g.delete(s)):g.delete(s),g.size===0)))break;d.recordUnknownSegments&&(a=oe.findPosition(u,t),a.marker=i,this.unknownSegments.push(a)),t+=r+1}else if(wa(i)){if(r=u.getUint16(t+2),i===218&&d.stopAfterSos!==!1)return;d.recordJpegSegments&&this.jpegSegments.push({offset:t,length:r,marker:i}),t+=r+1}}return t}mergeMultiSegments(){if(!this.appSegments.some(e=>e.multiSegment))return;let t=function(e,i){let r,s,o,a=new Map;for(let l=0;l<e.length;l++)r=e[l],s=r[i],a.has(s)?o=a.get(s):a.set(s,o=[]),o.push(r);return Array.from(a)}(this.appSegments,"type");this.mergedAppSegments=t.map(([e,i])=>{let r=Ct.get(e,this.options);return r.handleMultiSegments?{type:e,chunk:r.handleMultiSegments(i)}:i[0]})}getSegment(t){return this.appSegments.find(e=>e.type===t)}async getOrFindSegment(t){let e=this.getSegment(t);return e===void 0&&(await this.findAppSegments(0,[t]),e=this.getSegment(t)),e}}st(kr,"type","jpeg"),be.set("jpeg",kr);const Pa=[void 0,1,1,2,4,8,1,1,2,4,8,4,8,4];class Ca extends oe{parseHeader(){var t=this.chunk.getUint16();t===18761?this.le=!0:t===19789&&(this.le=!1),this.chunk.le=this.le,this.headerParsed=!0}parseTags(t,e,i=new Map){let{pick:r,skip:s}=this.options[e];r=new Set(r);let o=r.size>0,a=s.size===0,l=this.chunk.getUint16(t);t+=2;for(let u=0;u<l;u++){let h=this.chunk.getUint16(t);if(o){if(r.has(h)&&(i.set(h,this.parseTag(t,h,e)),r.delete(h),r.size===0))break}else!a&&s.has(h)||i.set(h,this.parseTag(t,h,e));t+=12}return i}parseTag(t,e,i){let{chunk:r}=this,s=r.getUint16(t+2),o=r.getUint32(t+4),a=Pa[s];if(a*o<=4?t+=8:t=r.getUint32(t+8),(s<1||s>13)&&kt(`Invalid TIFF value type. block: ${i.toUpperCase()}, tag: ${e.toString(16)}, type: ${s}, offset ${t}`),t>r.byteLength&&kt(`Invalid TIFF value offset. block: ${i.toUpperCase()}, tag: ${e.toString(16)}, type: ${s}, offset ${t} is outside of chunk size ${r.byteLength}`),s===1)return r.getUint8Array(t,o);if(s===2)return Ze(r.getString(t,o));if(s===7)return r.getUint8Array(t,o);if(o===1)return this.parseTagValue(s,t);{let l=new(function(h){switch(h){case 1:return Uint8Array;case 3:return Uint16Array;case 4:return Uint32Array;case 5:return Array;case 6:return Int8Array;case 8:return Int16Array;case 9:return Int32Array;case 10:return Array;case 11:return Float32Array;case 12:return Float64Array;default:return Array}}(s))(o),u=a;for(let h=0;h<o;h++)l[h]=this.parseTagValue(s,t),t+=u;return l}}parseTagValue(t,e){let{chunk:i}=this;switch(t){case 1:return i.getUint8(e);case 3:return i.getUint16(e);case 4:return i.getUint32(e);case 5:return i.getUint32(e)/i.getUint32(e+4);case 6:return i.getInt8(e);case 8:return i.getInt16(e);case 9:return i.getInt32(e);case 10:return i.getInt32(e)/i.getInt32(e+4);case 11:return i.getFloat(e);case 12:return i.getDouble(e);case 13:return i.getUint32(e);default:kt(`Invalid tiff type ${t}`)}}}class Ni extends Ca{static canHandle(t,e){return t.getUint8(e+1)===225&&t.getUint32(e+4)===1165519206&&t.getUint16(e+8)===0}async parse(){this.parseHeader();let{options:t}=this;return t.ifd0.enabled&&await this.parseIfd0Block(),t.exif.enabled&&await this.safeParse("parseExifBlock"),t.gps.enabled&&await this.safeParse("parseGpsBlock"),t.interop.enabled&&await this.safeParse("parseInteropBlock"),t.ifd1.enabled&&await this.safeParse("parseThumbnailBlock"),this.createOutput()}safeParse(t){let e=this[t]();return e.catch!==void 0&&(e=e.catch(this.handleError)),e}findIfd0Offset(){this.ifd0Offset===void 0&&(this.ifd0Offset=this.chunk.getUint32(4))}findIfd1Offset(){if(this.ifd1Offset===void 0){this.findIfd0Offset();let t=this.chunk.getUint16(this.ifd0Offset),e=this.ifd0Offset+2+12*t;this.ifd1Offset=this.chunk.getUint32(e)}}parseBlock(t,e){let i=new Map;return this[e]=i,this.parseTags(t,e,i),i}async parseIfd0Block(){if(this.ifd0)return;let{file:t}=this;this.findIfd0Offset(),this.ifd0Offset<8&&kt("Malformed EXIF data"),!t.chunked&&this.ifd0Offset>t.byteLength&&kt(`IFD0 offset points to outside of file.
this.ifd0Offset: ${this.ifd0Offset}, file.byteLength: ${t.byteLength}`),t.tiff&&await t.ensureChunk(this.ifd0Offset,Vi(this.options));let e=this.parseBlock(this.ifd0Offset,"ifd0");return e.size!==0?(this.exifOffset=e.get(34665),this.interopOffset=e.get(40965),this.gpsOffset=e.get(34853),this.xmp=e.get(700),this.iptc=e.get(33723),this.icc=e.get(34675),this.options.sanitize&&(e.delete(34665),e.delete(40965),e.delete(34853),e.delete(700),e.delete(33723),e.delete(34675)),e):void 0}async parseExifBlock(){if(this.exif||(this.ifd0||await this.parseIfd0Block(),this.exifOffset===void 0))return;this.file.tiff&&await this.file.ensureChunk(this.exifOffset,Vi(this.options));let t=this.parseBlock(this.exifOffset,"exif");return this.interopOffset||(this.interopOffset=t.get(40965)),this.makerNote=t.get(37500),this.userComment=t.get(37510),this.options.sanitize&&(t.delete(40965),t.delete(37500),t.delete(37510)),this.unpack(t,41728),this.unpack(t,41729),t}unpack(t,e){let i=t.get(e);i&&i.length===1&&t.set(e,i[0])}async parseGpsBlock(){if(this.gps||(this.ifd0||await this.parseIfd0Block(),this.gpsOffset===void 0))return;let t=this.parseBlock(this.gpsOffset,"gps");return t&&t.has(2)&&t.has(4)&&(t.set("latitude",Fr(...t.get(2),t.get(1))),t.set("longitude",Fr(...t.get(4),t.get(3)))),t}async parseInteropBlock(){if(!this.interop&&(this.ifd0||await this.parseIfd0Block(),this.interopOffset!==void 0||this.exif||await this.parseExifBlock(),this.interopOffset!==void 0))return this.parseBlock(this.interopOffset,"interop")}async parseThumbnailBlock(t=!1){if(!this.ifd1&&!this.ifd1Parsed&&(!this.options.mergeOutput||t))return this.findIfd1Offset(),this.ifd1Offset>0&&(this.parseBlock(this.ifd1Offset,"ifd1"),this.ifd1Parsed=!0),this.ifd1}async extractThumbnail(){if(this.headerParsed||this.parseHeader(),this.ifd1Parsed||await this.parseThumbnailBlock(!0),this.ifd1===void 0)return;let t=this.ifd1.get(513),e=this.ifd1.get(514);return this.chunk.getUint8Array(t,e)}get image(){return this.ifd0}get thumbnail(){return this.ifd1}createOutput(){let t,e,i,r={};for(e of St)if(t=this[e],!is(t))if(i=this.canTranslate?this.translateBlock(t,e):Object.fromEntries(t),this.options.mergeOutput){if(e==="ifd1")continue;Object.assign(r,i)}else r[e]=i;return this.makerNote&&(r.makerNote=this.makerNote),this.userComment&&(r.userComment=this.userComment),r}assignToOutput(t,e){if(this.globalOptions.mergeOutput)Object.assign(t,e);else for(let[i,r]of Object.entries(e))this.assignObjectToOutput(t,i,r)}}function Fr(n,t,e,i){var r=n+t/60+e/3600;return i!=="S"&&i!=="W"||(r*=-1),r}st(Ni,"type","tiff"),st(Ni,"headerLength",10),Ct.set("tiff",Ni);var ka=Object.freeze({__proto__:null,default:Ma,Exifr:De,fileParsers:be,segmentParsers:Ct,fileReaders:Me,tagKeys:Rt,tagValues:le,tagRevivers:Be,createDictionary:Pt,extendDictionary:cn,fetchUrlAsArrayBuffer:on,readBlobAsArrayBuffer:ln,chunkedProps:Re,otherSegments:bn,segments:un,tiffBlocks:St,segmentsAndBlocks:Ne,tiffExtractables:Le,inheritables:_n,allFormatters:Ee,Options:hn,parse:ui});const $i={ifd0:!1,ifd1:!1,exif:!1,gps:!1,interop:!1,sanitize:!1,reviveValues:!0,translateKeys:!1,translateValues:!1,mergeOutput:!1},Ji=Object.assign({},$i,{firstChunkSize:4e4,gps:[1,2,3,4]});async function as(n){let t=new De(Ji);await t.read(n);let e=await t.parse();if(e&&e.gps){let{latitude:i,longitude:r}=e.gps;return{latitude:i,longitude:r}}}const Zi=Object.assign({},$i,{tiff:!1,ifd1:!0,mergeOutput:!1});async function os(n){let t=new De(Zi);await t.read(n);let e=await t.extractThumbnail();return e&&ci?li.from(e):e}async function ls(n){let t=await this.thumbnail(n);if(t!==void 0){let e=new Blob([t]);return URL.createObjectURL(e)}}const tr=Object.assign({},$i,{firstChunkSize:4e4,ifd0:[274]});async function er(n){let t=new De(tr);await t.read(n);let e=await t.parse();if(e&&e.ifd0)return e.ifd0[274]}const nr=Object.freeze({1:{dimensionSwapped:!1,scaleX:1,scaleY:1,deg:0,rad:0},2:{dimensionSwapped:!1,scaleX:-1,scaleY:1,deg:0,rad:0},3:{dimensionSwapped:!1,scaleX:1,scaleY:1,deg:180,rad:180*Math.PI/180},4:{dimensionSwapped:!1,scaleX:-1,scaleY:1,deg:180,rad:180*Math.PI/180},5:{dimensionSwapped:!0,scaleX:1,scaleY:-1,deg:90,rad:90*Math.PI/180},6:{dimensionSwapped:!0,scaleX:1,scaleY:1,deg:90,rad:90*Math.PI/180},7:{dimensionSwapped:!0,scaleX:1,scaleY:-1,deg:270,rad:270*Math.PI/180},8:{dimensionSwapped:!0,scaleX:1,scaleY:1,deg:270,rad:270*Math.PI/180}});let tn=!0,en=!0;if(typeof navigator=="object"){let n=navigator.userAgent;if(n.includes("iPad")||n.includes("iPhone")){let t=n.match(/OS (\d+)_(\d+)/);if(t){let[,e,i]=t;tn=Number(e)+.1*Number(i)<13.4,en=!1}}else if(n.includes("OS X 10")){let[,t]=n.match(/OS X 10[_.](\d+)/);tn=en=Number(t)<15}if(n.includes("Chrome/")){let[,t]=n.match(/Chrome\/(\d+)/);tn=en=Number(t)<81}else if(n.includes("Firefox/")){let[,t]=n.match(/Firefox\/(\d+)/);tn=en=Number(t)<77}}async function cs(n){let t=await er(n);return Object.assign({canvas:tn,css:en},nr[t])}class Fa extends jt{constructor(...t){super(...t),st(this,"ranges",new Aa),this.byteLength!==0&&this.ranges.add(0,this.byteLength)}_tryExtend(t,e,i){if(t===0&&this.byteLength===0&&i){let r=new DataView(i.buffer||i,i.byteOffset,i.byteLength);this._swapDataView(r)}else{let r=t+e;if(r>this.byteLength){let{dataView:s}=this._extend(r);this._swapDataView(s)}}}_extend(t){let e;e=ci?li.allocUnsafe(t):new Uint8Array(t);let i=new DataView(e.buffer,e.byteOffset,e.byteLength);return e.set(new Uint8Array(this.buffer,this.byteOffset,this.byteLength),0),{uintView:e,dataView:i}}subarray(t,e,i=!1){return e=e||this._lengthToEnd(t),i&&this._tryExtend(t,e),this.ranges.add(t,e),super.subarray(t,e)}set(t,e,i=!1){i&&this._tryExtend(e,t.byteLength,t);let r=super.set(t,e);return this.ranges.add(e,r.byteLength),r}async ensureChunk(t,e){this.chunked&&(this.ranges.available(t,e)||await this.readChunk(t,e))}available(t,e){return this.ranges.available(t,e)}}class Aa{constructor(){st(this,"list",[])}get length(){return this.list.length}add(t,e,i=0){let r=t+e,s=this.list.filter(o=>Ar(t,o.offset,r)||Ar(t,o.end,r));if(s.length>0){t=Math.min(t,...s.map(a=>a.offset)),r=Math.max(r,...s.map(a=>a.end)),e=r-t;let o=s.shift();o.offset=t,o.length=e,o.end=r,this.list=this.list.filter(a=>!s.includes(a))}else this.list.push({offset:t,length:e,end:r})}available(t,e){let i=t+e;return this.list.some(r=>r.offset<=t&&i<=r.end)}}function Ar(n,t,e){return n<=t&&t<=e}class di extends Fa{constructor(t,e){super(0),st(this,"chunksRead",0),this.input=t,this.options=e}async readWhole(){this.chunked=!1,await this.readChunk(this.nextChunkOffset)}async readChunked(){this.chunked=!0,await this.readChunk(0,this.options.firstChunkSize)}async readNextChunk(t=this.nextChunkOffset){if(this.fullyRead)return this.chunksRead++,!1;let e=this.options.chunkSize,i=await this.readChunk(t,e);return!!i&&i.byteLength===e}async readChunk(t,e){if(this.chunksRead++,(e=this.safeWrapAddress(t,e))!==0)return this._readChunk(t,e)}safeWrapAddress(t,e){return this.size!==void 0&&t+e>this.size?Math.max(0,this.size-t):e}get nextChunkOffset(){if(this.ranges.list.length!==0)return this.ranges.list[0].length}get canReadNextChunk(){return this.chunksRead<this.options.chunkLimit}get fullyRead(){return this.size!==void 0&&this.nextChunkOffset===this.size}read(){return this.options.chunked?this.readChunked():this.readWhole()}close(){}}Me.set("blob",class extends di{async readWhole(){this.chunked=!1;let n=await ln(this.input);this._swapArrayBuffer(n)}readChunked(){return this.chunked=!0,this.size=this.input.size,super.readChunked()}async _readChunk(n,t){let e=t?n+t:void 0,i=this.input.slice(n,e),r=await ln(i);return this.set(r,n,!0)}});var Ta=Object.freeze({__proto__:null,default:ka,Exifr:De,fileParsers:be,segmentParsers:Ct,fileReaders:Me,tagKeys:Rt,tagValues:le,tagRevivers:Be,createDictionary:Pt,extendDictionary:cn,fetchUrlAsArrayBuffer:on,readBlobAsArrayBuffer:ln,chunkedProps:Re,otherSegments:bn,segments:un,tiffBlocks:St,segmentsAndBlocks:Ne,tiffExtractables:Le,inheritables:_n,allFormatters:Ee,Options:hn,parse:ui,gpsOnlyOptions:Ji,gps:as,thumbnailOnlyOptions:Zi,thumbnail:os,thumbnailUrl:ls,orientationOnlyOptions:tr,orientation:er,rotations:nr,get rotateCanvas(){return tn},get rotateCss(){return en},rotation:cs});Me.set("url",class extends di{async readWhole(){this.chunked=!1;let n=await on(this.input);n instanceof ArrayBuffer?this._swapArrayBuffer(n):n instanceof Uint8Array&&this._swapBuffer(n)}async _readChunk(n,t){let e=t?n+t-1:void 0,i=this.options.httpHeaders||{};(n||e)&&(i.range=`bytes=${[n,e].join("-")}`);let r=await qi(this.input,{headers:i}),s=await r.arrayBuffer(),o=s.byteLength;if(r.status!==416)return o!==t&&(this.size=n+o),this.set(s,n,!0)}});jt.prototype.getUint64=function(n){let t=this.getUint32(n),e=this.getUint32(n+4);return t<1048575?t<<32|e:typeof Yn!==void 0?(console.warn("Using BigInt because of type 64uint but JS can only handle 53b numbers."),Yn(t)<<Yn(32)|Yn(e)):void kt("Trying to read 64b value but JS can only handle 53b numbers.")};class Ia extends hi{parseBoxes(t=0){let e=[];for(;t<this.file.byteLength-4;){let i=this.parseBoxHead(t);if(e.push(i),i.length===0)break;t+=i.length}return e}parseSubBoxes(t){t.boxes=this.parseBoxes(t.start)}findBox(t,e){return t.boxes===void 0&&this.parseSubBoxes(t),t.boxes.find(i=>i.kind===e)}parseBoxHead(t){let e=this.file.getUint32(t),i=this.file.getString(t+4,4),r=t+8;return e===1&&(e=this.file.getUint64(t+8),r+=8),{offset:t,length:e,kind:i,start:r}}parseBoxFullHead(t){if(t.version!==void 0)return;let e=this.file.getUint32(t.start);t.version=e>>24,t.start+=4}}class us extends Ia{static canHandle(t,e){if(e!==0)return!1;let i=t.getUint16(2);if(i>50)return!1;let r=16,s=[];for(;r<i;)s.push(t.getString(r,4)),r+=4;return s.includes(this.type)}async parse(){let t=this.file.getUint32(0),e=this.parseBoxHead(t);for(;e.kind!=="meta";)t+=e.length,await this.file.ensureChunk(t,16),e=this.parseBoxHead(t);await this.file.ensureChunk(e.offset,e.length),this.parseBoxFullHead(e),this.parseSubBoxes(e),this.options.icc.enabled&&await this.findIcc(e),this.options.tiff.enabled&&await this.findExif(e)}async registerSegment(t,e,i){await this.file.ensureChunk(e,i);let r=this.file.subarray(e,i);this.createParser(t,r)}async findIcc(t){let e=this.findBox(t,"iprp");if(e===void 0)return;let i=this.findBox(e,"ipco");if(i===void 0)return;let r=this.findBox(i,"colr");r!==void 0&&await this.registerSegment("icc",r.offset+12,r.length)}async findExif(t){let e=this.findBox(t,"iinf");if(e===void 0)return;let i=this.findBox(t,"iloc");if(i===void 0)return;let r=this.findExifLocIdInIinf(e),s=this.findExtentInIloc(i,r);if(s===void 0)return;let[o,a]=s;await this.file.ensureChunk(o,a);let l=4+this.file.getUint32(o);o+=l,a-=l,await this.registerSegment("tiff",o,a)}findExifLocIdInIinf(t){this.parseBoxFullHead(t);let e,i,r,s,o=t.start,a=this.file.getUint16(o);for(o+=2;a--;){if(e=this.parseBoxHead(o),this.parseBoxFullHead(e),i=e.start,e.version>=2&&(r=e.version===3?4:2,s=this.file.getString(i+r+2,4),s==="Exif"))return this.file.getUintBytes(i,r);o+=e.length}}get8bits(t){let e=this.file.getUint8(t);return[e>>4,15&e]}findExtentInIloc(t,e){this.parseBoxFullHead(t);let i=t.start,[r,s]=this.get8bits(i++),[o,a]=this.get8bits(i++),l=t.version===2?4:2,u=t.version===1||t.version===2?2:0,h=a+r+s,c=t.version===2?4:2,g=this.file.getUintBytes(i,c);for(i+=c;g--;){let d=this.file.getUintBytes(i,l);i+=l+u+2+o;let f=this.file.getUint16(i);if(i+=2,d===e)return f>1&&console.warn(`ILOC box has more than one extent but we're only processing one
Please create an issue at https://github.com/MikeKovarik/exifr with this file`),[this.file.getUintBytes(i+a,r),this.file.getUintBytes(i+a+r,s)];i+=f*h}}}class hs extends us{}st(hs,"type","heic");class Tr extends us{}st(Tr,"type","avif"),be.set("heic",hs),be.set("avif",Tr),Pt(Rt,["ifd0","ifd1"],[[256,"ImageWidth"],[257,"ImageHeight"],[258,"BitsPerSample"],[259,"Compression"],[262,"PhotometricInterpretation"],[270,"ImageDescription"],[271,"Make"],[272,"Model"],[273,"StripOffsets"],[274,"Orientation"],[277,"SamplesPerPixel"],[278,"RowsPerStrip"],[279,"StripByteCounts"],[282,"XResolution"],[283,"YResolution"],[284,"PlanarConfiguration"],[296,"ResolutionUnit"],[301,"TransferFunction"],[305,"Software"],[306,"ModifyDate"],[315,"Artist"],[316,"HostComputer"],[317,"Predictor"],[318,"WhitePoint"],[319,"PrimaryChromaticities"],[513,"ThumbnailOffset"],[514,"ThumbnailLength"],[529,"YCbCrCoefficients"],[530,"YCbCrSubSampling"],[531,"YCbCrPositioning"],[532,"ReferenceBlackWhite"],[700,"ApplicationNotes"],[33432,"Copyright"],[33723,"IPTC"],[34665,"ExifIFD"],[34675,"ICC"],[34853,"GpsIFD"],[330,"SubIFD"],[40965,"InteropIFD"],[40091,"XPTitle"],[40092,"XPComment"],[40093,"XPAuthor"],[40094,"XPKeywords"],[40095,"XPSubject"]]),Pt(Rt,"exif",[[33434,"ExposureTime"],[33437,"FNumber"],[34850,"ExposureProgram"],[34852,"SpectralSensitivity"],[34855,"ISO"],[34858,"TimeZoneOffset"],[34859,"SelfTimerMode"],[34864,"SensitivityType"],[34865,"StandardOutputSensitivity"],[34866,"RecommendedExposureIndex"],[34867,"ISOSpeed"],[34868,"ISOSpeedLatitudeyyy"],[34869,"ISOSpeedLatitudezzz"],[36864,"ExifVersion"],[36867,"DateTimeOriginal"],[36868,"CreateDate"],[36873,"GooglePlusUploadCode"],[36880,"OffsetTime"],[36881,"OffsetTimeOriginal"],[36882,"OffsetTimeDigitized"],[37121,"ComponentsConfiguration"],[37122,"CompressedBitsPerPixel"],[37377,"ShutterSpeedValue"],[37378,"ApertureValue"],[37379,"BrightnessValue"],[37380,"ExposureCompensation"],[37381,"MaxApertureValue"],[37382,"SubjectDistance"],[37383,"MeteringMode"],[37384,"LightSource"],[37385,"Flash"],[37386,"FocalLength"],[37393,"ImageNumber"],[37394,"SecurityClassification"],[37395,"ImageHistory"],[37396,"SubjectArea"],[37500,"MakerNote"],[37510,"UserComment"],[37520,"SubSecTime"],[37521,"SubSecTimeOriginal"],[37522,"SubSecTimeDigitized"],[37888,"AmbientTemperature"],[37889,"Humidity"],[37890,"Pressure"],[37891,"WaterDepth"],[37892,"Acceleration"],[37893,"CameraElevationAngle"],[40960,"FlashpixVersion"],[40961,"ColorSpace"],[40962,"ExifImageWidth"],[40963,"ExifImageHeight"],[40964,"RelatedSoundFile"],[41483,"FlashEnergy"],[41486,"FocalPlaneXResolution"],[41487,"FocalPlaneYResolution"],[41488,"FocalPlaneResolutionUnit"],[41492,"SubjectLocation"],[41493,"ExposureIndex"],[41495,"SensingMethod"],[41728,"FileSource"],[41729,"SceneType"],[41730,"CFAPattern"],[41985,"CustomRendered"],[41986,"ExposureMode"],[41987,"WhiteBalance"],[41988,"DigitalZoomRatio"],[41989,"FocalLengthIn35mmFormat"],[41990,"SceneCaptureType"],[41991,"GainControl"],[41992,"Contrast"],[41993,"Saturation"],[41994,"Sharpness"],[41996,"SubjectDistanceRange"],[42016,"ImageUniqueID"],[42032,"OwnerName"],[42033,"SerialNumber"],[42034,"LensInfo"],[42035,"LensMake"],[42036,"LensModel"],[42037,"LensSerialNumber"],[42080,"CompositeImage"],[42081,"CompositeImageCount"],[42082,"CompositeImageExposureTimes"],[42240,"Gamma"],[59932,"Padding"],[59933,"OffsetSchema"],[65e3,"OwnerName"],[65001,"SerialNumber"],[65002,"Lens"],[65100,"RawFile"],[65101,"Converter"],[65102,"WhiteBalance"],[65105,"Exposure"],[65106,"Shadows"],[65107,"Brightness"],[65108,"Contrast"],[65109,"Saturation"],[65110,"Sharpness"],[65111,"Smoothness"],[65112,"MoireFilter"],[40965,"InteropIFD"]]),Pt(Rt,"gps",[[0,"GPSVersionID"],[1,"GPSLatitudeRef"],[2,"GPSLatitude"],[3,"GPSLongitudeRef"],[4,"GPSLongitude"],[5,"GPSAltitudeRef"],[6,"GPSAltitude"],[7,"GPSTimeStamp"],[8,"GPSSatellites"],[9,"GPSStatus"],[10,"GPSMeasureMode"],[11,"GPSDOP"],[12,"GPSSpeedRef"],[13,"GPSSpeed"],[14,"GPSTrackRef"],[15,"GPSTrack"],[16,"GPSImgDirectionRef"],[17,"GPSImgDirection"],[18,"GPSMapDatum"],[19,"GPSDestLatitudeRef"],[20,"GPSDestLatitude"],[21,"GPSDestLongitudeRef"],[22,"GPSDestLongitude"],[23,"GPSDestBearingRef"],[24,"GPSDestBearing"],[25,"GPSDestDistanceRef"],[26,"GPSDestDistance"],[27,"GPSProcessingMethod"],[28,"GPSAreaInformation"],[29,"GPSDateStamp"],[30,"GPSDifferential"],[31,"GPSHPositioningError"]]),Pt(le,["ifd0","ifd1"],[[274,{1:"Horizontal (normal)",2:"Mirror horizontal",3:"Rotate 180",4:"Mirror vertical",5:"Mirror horizontal and rotate 270 CW",6:"Rotate 90 CW",7:"Mirror horizontal and rotate 90 CW",8:"Rotate 270 CW"}],[296,{1:"None",2:"inches",3:"cm"}]]);let An=Pt(le,"exif",[[34850,{0:"Not defined",1:"Manual",2:"Normal program",3:"Aperture priority",4:"Shutter priority",5:"Creative program",6:"Action program",7:"Portrait mode",8:"Landscape mode"}],[37121,{0:"-",1:"Y",2:"Cb",3:"Cr",4:"R",5:"G",6:"B"}],[37383,{0:"Unknown",1:"Average",2:"CenterWeightedAverage",3:"Spot",4:"MultiSpot",5:"Pattern",6:"Partial",255:"Other"}],[37384,{0:"Unknown",1:"Daylight",2:"Fluorescent",3:"Tungsten (incandescent light)",4:"Flash",9:"Fine weather",10:"Cloudy weather",11:"Shade",12:"Daylight fluorescent (D 5700 - 7100K)",13:"Day white fluorescent (N 4600 - 5400K)",14:"Cool white fluorescent (W 3900 - 4500K)",15:"White fluorescent (WW 3200 - 3700K)",17:"Standard light A",18:"Standard light B",19:"Standard light C",20:"D55",21:"D65",22:"D75",23:"D50",24:"ISO studio tungsten",255:"Other"}],[37385,{0:"Flash did not fire",1:"Flash fired",5:"Strobe return light not detected",7:"Strobe return light detected",9:"Flash fired, compulsory flash mode",13:"Flash fired, compulsory flash mode, return light not detected",15:"Flash fired, compulsory flash mode, return light detected",16:"Flash did not fire, compulsory flash mode",24:"Flash did not fire, auto mode",25:"Flash fired, auto mode",29:"Flash fired, auto mode, return light not detected",31:"Flash fired, auto mode, return light detected",32:"No flash function",65:"Flash fired, red-eye reduction mode",69:"Flash fired, red-eye reduction mode, return light not detected",71:"Flash fired, red-eye reduction mode, return light detected",73:"Flash fired, compulsory flash mode, red-eye reduction mode",77:"Flash fired, compulsory flash mode, red-eye reduction mode, return light not detected",79:"Flash fired, compulsory flash mode, red-eye reduction mode, return light detected",89:"Flash fired, auto mode, red-eye reduction mode",93:"Flash fired, auto mode, return light not detected, red-eye reduction mode",95:"Flash fired, auto mode, return light detected, red-eye reduction mode"}],[41495,{1:"Not defined",2:"One-chip color area sensor",3:"Two-chip color area sensor",4:"Three-chip color area sensor",5:"Color sequential area sensor",7:"Trilinear sensor",8:"Color sequential linear sensor"}],[41728,{1:"Film Scanner",2:"Reflection Print Scanner",3:"Digital Camera"}],[41729,{1:"Directly photographed"}],[41985,{0:"Normal",1:"Custom",2:"HDR (no original saved)",3:"HDR (original saved)",4:"Original (for HDR)",6:"Panorama",7:"Portrait HDR",8:"Portrait"}],[41986,{0:"Auto",1:"Manual",2:"Auto bracket"}],[41987,{0:"Auto",1:"Manual"}],[41990,{0:"Standard",1:"Landscape",2:"Portrait",3:"Night",4:"Other"}],[41991,{0:"None",1:"Low gain up",2:"High gain up",3:"Low gain down",4:"High gain down"}],[41996,{0:"Unknown",1:"Macro",2:"Close",3:"Distant"}],[42080,{0:"Unknown",1:"Not a Composite Image",2:"General Composite Image",3:"Composite Image Captured While Shooting"}]]);const Ir={1:"No absolute unit of measurement",2:"Inch",3:"Centimeter"};An.set(37392,Ir),An.set(41488,Ir);const Li={0:"Normal",1:"Low",2:"High"};function Rr(n){return typeof n=="object"&&n.length!==void 0?n[0]:n}function Nr(n){let t=Array.from(n).slice(1);return t[1]>15&&(t=t.map(e=>String.fromCharCode(e))),t[2]!=="0"&&t[2]!==0||t.pop(),t.join(".")}function Ei(n){if(typeof n=="string"){var[t,e,i,r,s,o]=n.trim().split(/[-: ]/g).map(Number),a=new Date(t,e-1,i);return Number.isNaN(r)||Number.isNaN(s)||Number.isNaN(o)||(a.setHours(r),a.setMinutes(s),a.setSeconds(o)),Number.isNaN(+a)?n:a}}function kn(n){if(typeof n=="string")return n;let t=[];if(n[1]===0&&n[n.length-1]===0)for(let e=0;e<n.length;e+=2)t.push(Lr(n[e+1],n[e]));else for(let e=0;e<n.length;e+=2)t.push(Lr(n[e],n[e+1]));return Ze(String.fromCodePoint(...t))}function Lr(n,t){return n<<8|t}An.set(41992,Li),An.set(41993,Li),An.set(41994,Li),Pt(Be,["ifd0","ifd1"],[[50827,function(n){return typeof n!="string"?rs(n):n}],[306,Ei],[40091,kn],[40092,kn],[40093,kn],[40094,kn],[40095,kn]]),Pt(Be,"exif",[[40960,Nr],[36864,Nr],[36867,Ei],[36868,Ei],[40962,Rr],[40963,Rr]]),Pt(Be,"gps",[[0,n=>Array.from(n).join(".")],[7,n=>Array.from(n).join(":")]]);class Ui extends oe{static canHandle(t,e){return t.getUint8(e+1)===225&&t.getUint32(e+4)===1752462448&&t.getString(e+4,20)==="http://ns.adobe.com/"}static headerLength(t,e){return t.getString(e+4,34)==="http://ns.adobe.com/xmp/extension/"?79:33}static findPosition(t,e){let i=super.findPosition(t,e);return i.multiSegment=i.extended=i.headerLength===79,i.multiSegment?(i.chunkCount=t.getUint8(e+72),i.chunkNumber=t.getUint8(e+76),t.getUint8(e+77)!==0&&i.chunkNumber++):(i.chunkCount=1/0,i.chunkNumber=-1),i}static handleMultiSegments(t){return t.map(e=>e.chunk.getString()).join("")}normalizeInput(t){return typeof t=="string"?t:jt.from(t).getString()}parse(t=this.chunk){if(!this.localOptions.parse)return t;t=function(s){let o={},a={};for(let l of ms)o[l]=[],a[l]=0;return s.replace(Ea,(l,u,h)=>{if(u==="<"){let c=++a[h];return o[h].push(c),`${l}#${c}`}return`${l}#${o[h].pop()}`})}(t);let e=mn.findAll(t,"rdf","Description");e.length===0&&e.push(new mn("rdf","Description",void 0,t));let i,r={};for(let s of e)for(let o of s.properties)i=La(o.ns,r),ds(o,i);return function(s){let o;for(let a in s)o=s[a]=ei(s[a]),o===void 0&&delete s[a];return ei(s)}(r)}assignToOutput(t,e){if(this.localOptions.parse)for(let[i,r]of Object.entries(e))switch(i){case"tiff":this.assignObjectToOutput(t,"ifd0",r);break;case"exif":this.assignObjectToOutput(t,"exif",r);break;case"xmlns":break;default:this.assignObjectToOutput(t,i,r)}else t.xmp=e}}st(Ui,"type","xmp"),st(Ui,"multiSegment",!0),Ct.set("xmp",Ui);class ii{static findAll(t){return fs(t,/([a-zA-Z0-9-]+):([a-zA-Z0-9-]+)=("[^"]*"|'[^']*')/gm).map(ii.unpackMatch)}static unpackMatch(t){let e=t[1],i=t[2],r=t[3].slice(1,-1);return r=ps(r),new ii(e,i,r)}constructor(t,e,i){this.ns=t,this.name=e,this.value=i}serialize(){return this.value}}class mn{static findAll(t,e,i){if(e!==void 0||i!==void 0){e=e||"[\\w\\d-]+",i=i||"[\\w\\d-]+";var r=new RegExp(`<(${e}):(${i})(#\\d+)?((\\s+?[\\w\\d-:]+=("[^"]*"|'[^']*'))*\\s*)(\\/>|>([\\s\\S]*?)<\\/\\1:\\2\\3>)`,"gm")}else r=/<([\w\d-]+):([\w\d-]+)(#\d+)?((\s+?[\w\d-:]+=("[^"]*"|'[^']*'))*\s*)(\/>|>([\s\S]*?)<\/\1:\2\3>)/gm;return fs(t,r).map(mn.unpackMatch)}static unpackMatch(t){let e=t[1],i=t[2],r=t[4],s=t[8];return new mn(e,i,r,s)}constructor(t,e,i,r){this.ns=t,this.name=e,this.attrString=i,this.innerXml=r,this.attrs=ii.findAll(i),this.children=mn.findAll(r),this.value=this.children.length===0?ps(r):void 0,this.properties=[...this.attrs,...this.children]}get isPrimitive(){return this.value!==void 0&&this.attrs.length===0&&this.children.length===0}get isListContainer(){return this.children.length===1&&this.children[0].isList}get isList(){let{ns:t,name:e}=this;return t==="rdf"&&(e==="Seq"||e==="Bag"||e==="Alt")}get isListItem(){return this.ns==="rdf"&&this.name==="li"}serialize(){if(this.properties.length===0&&this.value===void 0)return;if(this.isPrimitive)return this.value;if(this.isListContainer)return this.children[0].serialize();if(this.isList)return Na(this.children.map(Ra));if(this.isListItem&&this.children.length===1&&this.attrs.length===0)return this.children[0].serialize();let t={};for(let e of this.properties)ds(e,t);return this.value!==void 0&&(t.value=this.value),ei(t)}}function ds(n,t){let e=n.serialize();e!==void 0&&(t[n.name]=e)}var Ra=n=>n.serialize(),Na=n=>n.length===1?n[0]:n,La=(n,t)=>t[n]?t[n]:t[n]={};function fs(n,t){let e,i=[];if(!n)return i;for(;(e=t.exec(n))!==null;)i.push(e);return i}function ps(n){if(function(i){return i==null||i==="null"||i==="undefined"||i===""||i.trim()===""}(n))return;let t=Number(n);if(!Number.isNaN(t))return t;let e=n.toLowerCase();return e==="true"||e!=="false"&&n.trim()}const ms=["rdf:li","rdf:Seq","rdf:Bag","rdf:Alt","rdf:Description"],Ea=new RegExp(`(<|\\/)(${ms.join("|")})`,"g");var gs=Object.freeze({__proto__:null,default:Ta,Exifr:De,fileParsers:be,segmentParsers:Ct,fileReaders:Me,tagKeys:Rt,tagValues:le,tagRevivers:Be,createDictionary:Pt,extendDictionary:cn,fetchUrlAsArrayBuffer:on,readBlobAsArrayBuffer:ln,chunkedProps:Re,otherSegments:bn,segments:un,tiffBlocks:St,segmentsAndBlocks:Ne,tiffExtractables:Le,inheritables:_n,allFormatters:Ee,Options:hn,parse:ui,gpsOnlyOptions:Ji,gps:as,thumbnailOnlyOptions:Zi,thumbnail:os,thumbnailUrl:ls,orientationOnlyOptions:tr,orientation:er,rotations:nr,get rotateCanvas(){return tn},get rotateCss(){return en},rotation:cs});let Er=ti("fs",n=>n.promises);Me.set("fs",class extends di{async readWhole(){this.chunked=!1,this.fs=await Er;let n=await this.fs.readFile(this.input);this._swapBuffer(n)}async readChunked(){this.chunked=!0,this.fs=await Er,await this.open(),await this.readChunk(0,this.options.firstChunkSize)}async open(){this.fh===void 0&&(this.fh=await this.fs.open(this.input,"r"),this.size=(await this.fh.stat(this.input)).size)}async _readChunk(n,t){this.fh===void 0&&await this.open(),n+t>this.size&&(t=this.size-n);var e=this.subarray(n,t,!0);return await this.fh.read(e.dataView,0,t,n),e}async close(){if(this.fh){let n=this.fh;this.fh=void 0,await n.close()}}});Me.set("base64",class extends di{constructor(...n){super(...n),this.input=this.input.replace(/^data:([^;]+);base64,/gim,""),this.size=this.input.length/4*3,this.input.endsWith("==")?this.size-=2:this.input.endsWith("=")&&(this.size-=1)}async _readChunk(n,t){let e,i,r=this.input;n===void 0?(n=0,e=0,i=0):(e=4*Math.floor(n/3),i=n-e/4*3),t===void 0&&(t=this.size);let s=n+t,o=e+4*Math.ceil(s/3);r=r.slice(e,o);let a=Math.min(t,this.size-n);if(ci){let l=li.from(r,"base64").slice(i,i+a);return this.set(l,n,!0)}{let l=this.subarray(n,a,!0),u=atob(r),h=l.toUint8();for(let c=0;c<a;c++)h[c]=u.charCodeAt(i+c);return l}}});class Ur extends hi{static canHandle(t,e){return e===18761||e===19789}extendOptions(t){let{ifd0:e,xmp:i,iptc:r,icc:s}=t;i.enabled&&e.deps.add(700),r.enabled&&e.deps.add(33723),s.enabled&&e.deps.add(34675),e.finalizeFilters()}async parse(){let{tiff:t,xmp:e,iptc:i,icc:r}=this.options;if(t.enabled||e.enabled||i.enabled||r.enabled){let s=Math.max(Vi(this.options),this.options.chunkSize);await this.file.ensureChunk(0,s),this.createParser("tiff",this.file),this.parsers.tiff.parseHeader(),await this.parsers.tiff.parseIfd0Block(),this.adaptTiffPropAsSegment("xmp"),this.adaptTiffPropAsSegment("iptc"),this.adaptTiffPropAsSegment("icc")}}adaptTiffPropAsSegment(t){if(this.parsers.tiff[t]){let e=this.parsers.tiff[t];this.injectSegment(t,e)}}}st(Ur,"type","tiff"),be.set("tiff",Ur);let Ua=ti("zlib");const Ba=["ihdr","iccp","text","itxt","exif"];class Br extends hi{constructor(...t){super(...t),st(this,"catchError",e=>this.errors.push(e)),st(this,"metaChunks",[]),st(this,"unknownChunks",[])}static canHandle(t,e){return e===35152&&t.getUint32(0)===2303741511&&t.getUint32(4)===218765834}async parse(){let{file:t}=this;await this.findPngChunksInRange(8,t.byteLength),await this.readSegments(this.metaChunks),this.findIhdr(),this.parseTextChunks(),await this.findExif().catch(this.catchError),await this.findXmp().catch(this.catchError),await this.findIcc().catch(this.catchError)}async findPngChunksInRange(t,e){let{file:i}=this;for(;t<e;){let r=i.getUint32(t),s=i.getUint32(t+4),o=i.getString(t+4,4).toLowerCase(),a=r+4+4+4,l={type:o,offset:t,length:a,start:t+4+4,size:r,marker:s};Ba.includes(o)?this.metaChunks.push(l):this.unknownChunks.push(l),t+=a}}parseTextChunks(){let t=this.metaChunks.filter(e=>e.type==="text");for(let e of t){let[i,r]=this.file.getString(e.start,e.size).split("\0");this.injectKeyValToIhdr(i,r)}}injectKeyValToIhdr(t,e){let i=this.parsers.ihdr;i&&i.raw.set(t,e)}findIhdr(){let t=this.metaChunks.find(e=>e.type==="ihdr");t&&this.options.ihdr.enabled!==!1&&this.createParser("ihdr",t.chunk)}async findExif(){let t=this.metaChunks.find(e=>e.type==="exif");t&&this.injectSegment("tiff",t.chunk)}async findXmp(){let t=this.metaChunks.filter(e=>e.type==="itxt");for(let e of t)e.chunk.getString(0,17)==="XML:com.adobe.xmp"&&this.injectSegment("xmp",e.chunk)}async findIcc(){let t=this.metaChunks.find(a=>a.type==="iccp");if(!t)return;let{chunk:e}=t,i=e.getUint8Array(0,81),r=0;for(;r<80&&i[r]!==0;)r++;let s=r+2,o=e.getString(0,r);if(this.injectKeyValToIhdr("ProfileName",o),Zn){let a=await Ua,l=e.getUint8Array(s);l=a.inflateSync(l),this.injectSegment("icc",l)}}}st(Br,"type","png"),be.set("png",Br),Pt(Rt,"interop",[[1,"InteropIndex"],[2,"InteropVersion"],[4096,"RelatedImageFileFormat"],[4097,"RelatedImageWidth"],[4098,"RelatedImageHeight"]]),cn(Rt,"ifd0",[[11,"ProcessingSoftware"],[254,"SubfileType"],[255,"OldSubfileType"],[263,"Thresholding"],[264,"CellWidth"],[265,"CellLength"],[266,"FillOrder"],[269,"DocumentName"],[280,"MinSampleValue"],[281,"MaxSampleValue"],[285,"PageName"],[286,"XPosition"],[287,"YPosition"],[290,"GrayResponseUnit"],[297,"PageNumber"],[321,"HalftoneHints"],[322,"TileWidth"],[323,"TileLength"],[332,"InkSet"],[337,"TargetPrinter"],[18246,"Rating"],[18249,"RatingPercent"],[33550,"PixelScale"],[34264,"ModelTransform"],[34377,"PhotoshopSettings"],[50706,"DNGVersion"],[50707,"DNGBackwardVersion"],[50708,"UniqueCameraModel"],[50709,"LocalizedCameraModel"],[50736,"DNGLensInfo"],[50739,"ShadowScale"],[50740,"DNGPrivateData"],[33920,"IntergraphMatrix"],[33922,"ModelTiePoint"],[34118,"SEMInfo"],[34735,"GeoTiffDirectory"],[34736,"GeoTiffDoubleParams"],[34737,"GeoTiffAsciiParams"],[50341,"PrintIM"],[50721,"ColorMatrix1"],[50722,"ColorMatrix2"],[50723,"CameraCalibration1"],[50724,"CameraCalibration2"],[50725,"ReductionMatrix1"],[50726,"ReductionMatrix2"],[50727,"AnalogBalance"],[50728,"AsShotNeutral"],[50729,"AsShotWhiteXY"],[50730,"BaselineExposure"],[50731,"BaselineNoise"],[50732,"BaselineSharpness"],[50734,"LinearResponseLimit"],[50735,"CameraSerialNumber"],[50741,"MakerNoteSafety"],[50778,"CalibrationIlluminant1"],[50779,"CalibrationIlluminant2"],[50781,"RawDataUniqueID"],[50827,"OriginalRawFileName"],[50828,"OriginalRawFileData"],[50831,"AsShotICCProfile"],[50832,"AsShotPreProfileMatrix"],[50833,"CurrentICCProfile"],[50834,"CurrentPreProfileMatrix"],[50879,"ColorimetricReference"],[50885,"SRawType"],[50898,"PanasonicTitle"],[50899,"PanasonicTitle2"],[50931,"CameraCalibrationSig"],[50932,"ProfileCalibrationSig"],[50933,"ProfileIFD"],[50934,"AsShotProfileName"],[50936,"ProfileName"],[50937,"ProfileHueSatMapDims"],[50938,"ProfileHueSatMapData1"],[50939,"ProfileHueSatMapData2"],[50940,"ProfileToneCurve"],[50941,"ProfileEmbedPolicy"],[50942,"ProfileCopyright"],[50964,"ForwardMatrix1"],[50965,"ForwardMatrix2"],[50966,"PreviewApplicationName"],[50967,"PreviewApplicationVersion"],[50968,"PreviewSettingsName"],[50969,"PreviewSettingsDigest"],[50970,"PreviewColorSpace"],[50971,"PreviewDateTime"],[50972,"RawImageDigest"],[50973,"OriginalRawFileDigest"],[50981,"ProfileLookTableDims"],[50982,"ProfileLookTableData"],[51043,"TimeCodes"],[51044,"FrameRate"],[51058,"TStop"],[51081,"ReelName"],[51089,"OriginalDefaultFinalSize"],[51090,"OriginalBestQualitySize"],[51091,"OriginalDefaultCropSize"],[51105,"CameraLabel"],[51107,"ProfileHueSatMapEncoding"],[51108,"ProfileLookTableEncoding"],[51109,"BaselineExposureOffset"],[51110,"DefaultBlackRender"],[51111,"NewRawImageDigest"],[51112,"RawToPreviewGain"]]);let Dr=[[273,"StripOffsets"],[279,"StripByteCounts"],[288,"FreeOffsets"],[289,"FreeByteCounts"],[291,"GrayResponseCurve"],[292,"T4Options"],[293,"T6Options"],[300,"ColorResponseUnit"],[320,"ColorMap"],[324,"TileOffsets"],[325,"TileByteCounts"],[326,"BadFaxLines"],[327,"CleanFaxData"],[328,"ConsecutiveBadFaxLines"],[330,"SubIFD"],[333,"InkNames"],[334,"NumberofInks"],[336,"DotRange"],[338,"ExtraSamples"],[339,"SampleFormat"],[340,"SMinSampleValue"],[341,"SMaxSampleValue"],[342,"TransferRange"],[343,"ClipPath"],[344,"XClipPathUnits"],[345,"YClipPathUnits"],[346,"Indexed"],[347,"JPEGTables"],[351,"OPIProxy"],[400,"GlobalParametersIFD"],[401,"ProfileType"],[402,"FaxProfile"],[403,"CodingMethods"],[404,"VersionYear"],[405,"ModeNumber"],[433,"Decode"],[434,"DefaultImageColor"],[435,"T82Options"],[437,"JPEGTables"],[512,"JPEGProc"],[515,"JPEGRestartInterval"],[517,"JPEGLosslessPredictors"],[518,"JPEGPointTransforms"],[519,"JPEGQTables"],[520,"JPEGDCTables"],[521,"JPEGACTables"],[559,"StripRowCounts"],[999,"USPTOMiscellaneous"],[18247,"XP_DIP_XML"],[18248,"StitchInfo"],[28672,"SonyRawFileType"],[28688,"SonyToneCurve"],[28721,"VignettingCorrection"],[28722,"VignettingCorrParams"],[28724,"ChromaticAberrationCorrection"],[28725,"ChromaticAberrationCorrParams"],[28726,"DistortionCorrection"],[28727,"DistortionCorrParams"],[29895,"SonyCropTopLeft"],[29896,"SonyCropSize"],[32781,"ImageID"],[32931,"WangTag1"],[32932,"WangAnnotation"],[32933,"WangTag3"],[32934,"WangTag4"],[32953,"ImageReferencePoints"],[32954,"RegionXformTackPoint"],[32955,"WarpQuadrilateral"],[32956,"AffineTransformMat"],[32995,"Matteing"],[32996,"DataType"],[32997,"ImageDepth"],[32998,"TileDepth"],[33300,"ImageFullWidth"],[33301,"ImageFullHeight"],[33302,"TextureFormat"],[33303,"WrapModes"],[33304,"FovCot"],[33305,"MatrixWorldToScreen"],[33306,"MatrixWorldToCamera"],[33405,"Model2"],[33421,"CFARepeatPatternDim"],[33422,"CFAPattern2"],[33423,"BatteryLevel"],[33424,"KodakIFD"],[33445,"MDFileTag"],[33446,"MDScalePixel"],[33447,"MDColorTable"],[33448,"MDLabName"],[33449,"MDSampleInfo"],[33450,"MDPrepDate"],[33451,"MDPrepTime"],[33452,"MDFileUnits"],[33589,"AdventScale"],[33590,"AdventRevision"],[33628,"UIC1Tag"],[33629,"UIC2Tag"],[33630,"UIC3Tag"],[33631,"UIC4Tag"],[33918,"IntergraphPacketData"],[33919,"IntergraphFlagRegisters"],[33921,"INGRReserved"],[34016,"Site"],[34017,"ColorSequence"],[34018,"IT8Header"],[34019,"RasterPadding"],[34020,"BitsPerRunLength"],[34021,"BitsPerExtendedRunLength"],[34022,"ColorTable"],[34023,"ImageColorIndicator"],[34024,"BackgroundColorIndicator"],[34025,"ImageColorValue"],[34026,"BackgroundColorValue"],[34027,"PixelIntensityRange"],[34028,"TransparencyIndicator"],[34029,"ColorCharacterization"],[34030,"HCUsage"],[34031,"TrapIndicator"],[34032,"CMYKEquivalent"],[34152,"AFCP_IPTC"],[34232,"PixelMagicJBIGOptions"],[34263,"JPLCartoIFD"],[34306,"WB_GRGBLevels"],[34310,"LeafData"],[34687,"TIFF_FXExtensions"],[34688,"MultiProfiles"],[34689,"SharedData"],[34690,"T88Options"],[34732,"ImageLayer"],[34750,"JBIGOptions"],[34856,"Opto-ElectricConvFactor"],[34857,"Interlace"],[34908,"FaxRecvParams"],[34909,"FaxSubAddress"],[34910,"FaxRecvTime"],[34929,"FedexEDR"],[34954,"LeafSubIFD"],[37387,"FlashEnergy"],[37388,"SpatialFrequencyResponse"],[37389,"Noise"],[37390,"FocalPlaneXResolution"],[37391,"FocalPlaneYResolution"],[37392,"FocalPlaneResolutionUnit"],[37397,"ExposureIndex"],[37398,"TIFF-EPStandardID"],[37399,"SensingMethod"],[37434,"CIP3DataFile"],[37435,"CIP3Sheet"],[37436,"CIP3Side"],[37439,"StoNits"],[37679,"MSDocumentText"],[37680,"MSPropertySetStorage"],[37681,"MSDocumentTextPosition"],[37724,"ImageSourceData"],[40965,"InteropIFD"],[40976,"SamsungRawPointersOffset"],[40977,"SamsungRawPointersLength"],[41217,"SamsungRawByteOrder"],[41218,"SamsungRawUnknown"],[41484,"SpatialFrequencyResponse"],[41485,"Noise"],[41489,"ImageNumber"],[41490,"SecurityClassification"],[41491,"ImageHistory"],[41494,"TIFF-EPStandardID"],[41995,"DeviceSettingDescription"],[42112,"GDALMetadata"],[42113,"GDALNoData"],[44992,"ExpandSoftware"],[44993,"ExpandLens"],[44994,"ExpandFilm"],[44995,"ExpandFilterLens"],[44996,"ExpandScanner"],[44997,"ExpandFlashLamp"],[46275,"HasselbladRawImage"],[48129,"PixelFormat"],[48130,"Transformation"],[48131,"Uncompressed"],[48132,"ImageType"],[48256,"ImageWidth"],[48257,"ImageHeight"],[48258,"WidthResolution"],[48259,"HeightResolution"],[48320,"ImageOffset"],[48321,"ImageByteCount"],[48322,"AlphaOffset"],[48323,"AlphaByteCount"],[48324,"ImageDataDiscard"],[48325,"AlphaDataDiscard"],[50215,"OceScanjobDesc"],[50216,"OceApplicationSelector"],[50217,"OceIDNumber"],[50218,"OceImageLogic"],[50255,"Annotations"],[50459,"HasselbladExif"],[50547,"OriginalFileName"],[50560,"USPTOOriginalContentType"],[50656,"CR2CFAPattern"],[50710,"CFAPlaneColor"],[50711,"CFALayout"],[50712,"LinearizationTable"],[50713,"BlackLevelRepeatDim"],[50714,"BlackLevel"],[50715,"BlackLevelDeltaH"],[50716,"BlackLevelDeltaV"],[50717,"WhiteLevel"],[50718,"DefaultScale"],[50719,"DefaultCropOrigin"],[50720,"DefaultCropSize"],[50733,"BayerGreenSplit"],[50737,"ChromaBlurRadius"],[50738,"AntiAliasStrength"],[50752,"RawImageSegmentation"],[50780,"BestQualityScale"],[50784,"AliasLayerMetadata"],[50829,"ActiveArea"],[50830,"MaskedAreas"],[50935,"NoiseReductionApplied"],[50974,"SubTileBlockSize"],[50975,"RowInterleaveFactor"],[51008,"OpcodeList1"],[51009,"OpcodeList2"],[51022,"OpcodeList3"],[51041,"NoiseProfile"],[51114,"CacheVersion"],[51125,"DefaultUserCrop"],[51157,"NikonNEFInfo"],[65024,"KdcIFD"]];cn(Rt,"ifd0",Dr),cn(Rt,"exif",Dr),Pt(le,"gps",[[23,{M:"Magnetic North",T:"True North"}],[25,{K:"Kilometers",M:"Miles",N:"Nautical Miles"}]]);class Bi extends oe{static canHandle(t,e){return t.getUint8(e+1)===224&&t.getUint32(e+4)===1246120262&&t.getUint8(e+8)===0}parse(){return this.parseTags(),this.translate(),this.output}parseTags(){this.raw=new Map([[0,this.chunk.getUint16(0)],[2,this.chunk.getUint8(2)],[3,this.chunk.getUint16(3)],[5,this.chunk.getUint16(5)],[7,this.chunk.getUint8(7)],[8,this.chunk.getUint8(8)]])}}st(Bi,"type","jfif"),st(Bi,"headerLength",9),Ct.set("jfif",Bi),Pt(Rt,"jfif",[[0,"JFIFVersion"],[2,"ResolutionUnit"],[3,"XResolution"],[5,"YResolution"],[7,"ThumbnailWidth"],[8,"ThumbnailHeight"]]);class Or extends oe{parse(){return this.parseTags(),this.translate(),this.output}parseTags(){this.raw=new Map([[0,this.chunk.getUint32(0)],[4,this.chunk.getUint32(4)],[8,this.chunk.getUint8(8)],[9,this.chunk.getUint8(9)],[10,this.chunk.getUint8(10)],[11,this.chunk.getUint8(11)],[12,this.chunk.getUint8(12)],...Array.from(this.raw)])}}st(Or,"type","ihdr"),Ct.set("ihdr",Or),Pt(Rt,"ihdr",[[0,"ImageWidth"],[4,"ImageHeight"],[8,"BitDepth"],[9,"ColorType"],[10,"Compression"],[11,"Filter"],[12,"Interlace"]]),Pt(le,"ihdr",[[9,{0:"Grayscale",2:"RGB",3:"Palette",4:"Grayscale with Alpha",6:"RGB with Alpha",DEFAULT:"Unknown"}],[10,{0:"Deflate/Inflate",DEFAULT:"Unknown"}],[11,{0:"Adaptive",DEFAULT:"Unknown"}],[12,{0:"Noninterlaced",1:"Adam7 Interlace",DEFAULT:"Unknown"}]]);class Jn extends oe{static canHandle(t,e){return t.getUint8(e+1)===226&&t.getUint32(e+4)===1229144927}static findPosition(t,e){let i=super.findPosition(t,e);return i.chunkNumber=t.getUint8(e+16),i.chunkCount=t.getUint8(e+17),i.multiSegment=i.chunkCount>1,i}static handleMultiSegments(t){return function(e){let i=function(r){let s=r[0].constructor,o=0;for(let u of r)o+=u.length;let a=new s(o),l=0;for(let u of r)a.set(u,l),l+=u.length;return a}(e.map(r=>r.chunk.toUint8()));return new jt(i)}(t)}parse(){return this.raw=new Map,this.parseHeader(),this.parseTags(),this.translate(),this.output}parseHeader(){let{raw:t}=this;this.chunk.byteLength<84&&kt("ICC header is too short");for(let[e,i]of Object.entries(Da)){e=parseInt(e,10);let r=i(this.chunk,e);r!=="\0\0\0\0"&&t.set(e,r)}}parseTags(){let t,e,i,r,s,{raw:o}=this,a=this.chunk.getUint32(128),l=132,u=this.chunk.byteLength;for(;a--;){if(t=this.chunk.getString(l,4),e=this.chunk.getUint32(l+4),i=this.chunk.getUint32(l+8),r=this.chunk.getString(e,4),e+i>u)return void console.warn("reached the end of the first ICC chunk. Enable options.tiff.multiSegment to read all ICC segments.");s=this.parseTag(r,e,i),s!==void 0&&s!=="\0\0\0\0"&&o.set(t,s),l+=12}}parseTag(t,e,i){switch(t){case"desc":return this.parseDesc(e);case"mluc":return this.parseMluc(e);case"text":return this.parseText(e,i);case"sig ":return this.parseSig(e)}if(!(e+i>this.chunk.byteLength))return this.chunk.getUint8Array(e,i)}parseDesc(t){let e=this.chunk.getUint32(t+8)-1;return Ze(this.chunk.getString(t+12,e))}parseText(t,e){return Ze(this.chunk.getString(t+8,e-8))}parseSig(t){return Ze(this.chunk.getString(t+8,4))}parseMluc(t){let{chunk:e}=this,i=e.getUint32(t+8),r=e.getUint32(t+12),s=t+16,o=[];for(let a=0;a<i;a++){let l=e.getString(s+0,2),u=e.getString(s+2,2),h=e.getUint32(s+4),c=e.getUint32(s+8)+t,g=Ze(e.getUnicodeString(c,h));o.push({lang:l,country:u,text:g}),s+=r}return i===1?o[0].text:o}translateValue(t,e){return typeof t=="string"?e[t]||e[t.toLowerCase()]||t:e[t]||t}}st(Jn,"type","icc"),st(Jn,"multiSegment",!0),st(Jn,"headerLength",18);const Da={4:ve,8:function(n,t){return[n.getUint8(t),n.getUint8(t+1)>>4,n.getUint8(t+1)%16].map(e=>e.toString(10)).join(".")},12:ve,16:ve,20:ve,24:function(n,t){const e=n.getUint16(t),i=n.getUint16(t+2)-1,r=n.getUint16(t+4),s=n.getUint16(t+6),o=n.getUint16(t+8),a=n.getUint16(t+10);return new Date(Date.UTC(e,i,r,s,o,a))},36:ve,40:ve,48:ve,52:ve,64:(n,t)=>n.getUint32(t),80:ve};function ve(n,t){return Ze(n.getString(t,4))}Ct.set("icc",Jn),Pt(Rt,"icc",[[4,"ProfileCMMType"],[8,"ProfileVersion"],[12,"ProfileClass"],[16,"ColorSpaceData"],[20,"ProfileConnectionSpace"],[24,"ProfileDateTime"],[36,"ProfileFileSignature"],[40,"PrimaryPlatform"],[44,"CMMFlags"],[48,"DeviceManufacturer"],[52,"DeviceModel"],[56,"DeviceAttributes"],[64,"RenderingIntent"],[68,"ConnectionSpaceIlluminant"],[80,"ProfileCreator"],[84,"ProfileID"],["Header","ProfileHeader"],["MS00","WCSProfiles"],["bTRC","BlueTRC"],["bXYZ","BlueMatrixColumn"],["bfd","UCRBG"],["bkpt","MediaBlackPoint"],["calt","CalibrationDateTime"],["chad","ChromaticAdaptation"],["chrm","Chromaticity"],["ciis","ColorimetricIntentImageState"],["clot","ColorantTableOut"],["clro","ColorantOrder"],["clrt","ColorantTable"],["cprt","ProfileCopyright"],["crdi","CRDInfo"],["desc","ProfileDescription"],["devs","DeviceSettings"],["dmdd","DeviceModelDesc"],["dmnd","DeviceMfgDesc"],["dscm","ProfileDescriptionML"],["fpce","FocalPlaneColorimetryEstimates"],["gTRC","GreenTRC"],["gXYZ","GreenMatrixColumn"],["gamt","Gamut"],["kTRC","GrayTRC"],["lumi","Luminance"],["meas","Measurement"],["meta","Metadata"],["mmod","MakeAndModel"],["ncl2","NamedColor2"],["ncol","NamedColor"],["ndin","NativeDisplayInfo"],["pre0","Preview0"],["pre1","Preview1"],["pre2","Preview2"],["ps2i","PS2RenderingIntent"],["ps2s","PostScript2CSA"],["psd0","PostScript2CRD0"],["psd1","PostScript2CRD1"],["psd2","PostScript2CRD2"],["psd3","PostScript2CRD3"],["pseq","ProfileSequenceDesc"],["psid","ProfileSequenceIdentifier"],["psvm","PS2CRDVMSize"],["rTRC","RedTRC"],["rXYZ","RedMatrixColumn"],["resp","OutputResponse"],["rhoc","ReflectionHardcopyOrigColorimetry"],["rig0","PerceptualRenderingIntentGamut"],["rig2","SaturationRenderingIntentGamut"],["rpoc","ReflectionPrintOutputColorimetry"],["sape","SceneAppearanceEstimates"],["scoe","SceneColorimetryEstimates"],["scrd","ScreeningDesc"],["scrn","Screening"],["targ","CharTarget"],["tech","Technology"],["vcgt","VideoCardGamma"],["view","ViewingConditions"],["vued","ViewingCondDesc"],["wtpt","MediaWhitePoint"]]);const Wn={"4d2p":"Erdt Systems",AAMA:"Aamazing Technologies",ACER:"Acer",ACLT:"Acolyte Color Research",ACTI:"Actix Sytems",ADAR:"Adara Technology",ADBE:"Adobe",ADI:"ADI Systems",AGFA:"Agfa Graphics",ALMD:"Alps Electric",ALPS:"Alps Electric",ALWN:"Alwan Color Expertise",AMTI:"Amiable Technologies",AOC:"AOC International",APAG:"Apago",APPL:"Apple Computer",AST:"AST","AT&T":"AT&T",BAEL:"BARBIERI electronic",BRCO:"Barco NV",BRKP:"Breakpoint",BROT:"Brother",BULL:"Bull",BUS:"Bus Computer Systems","C-IT":"C-Itoh",CAMR:"Intel",CANO:"Canon",CARR:"Carroll Touch",CASI:"Casio",CBUS:"Colorbus PL",CEL:"Crossfield",CELx:"Crossfield",CGS:"CGS Publishing Technologies International",CHM:"Rochester Robotics",CIGL:"Colour Imaging Group, London",CITI:"Citizen",CL00:"Candela",CLIQ:"Color IQ",CMCO:"Chromaco",CMiX:"CHROMiX",COLO:"Colorgraphic Communications",COMP:"Compaq",COMp:"Compeq/Focus Technology",CONR:"Conrac Display Products",CORD:"Cordata Technologies",CPQ:"Compaq",CPRO:"ColorPro",CRN:"Cornerstone",CTX:"CTX International",CVIS:"ColorVision",CWC:"Fujitsu Laboratories",DARI:"Darius Technology",DATA:"Dataproducts",DCP:"Dry Creek Photo",DCRC:"Digital Contents Resource Center, Chung-Ang University",DELL:"Dell Computer",DIC:"Dainippon Ink and Chemicals",DICO:"Diconix",DIGI:"Digital","DL&C":"Digital Light & Color",DPLG:"Doppelganger",DS:"Dainippon Screen",DSOL:"DOOSOL",DUPN:"DuPont",EPSO:"Epson",ESKO:"Esko-Graphics",ETRI:"Electronics and Telecommunications Research Institute",EVER:"Everex Systems",EXAC:"ExactCODE",Eizo:"Eizo",FALC:"Falco Data Products",FF:"Fuji Photo Film",FFEI:"FujiFilm Electronic Imaging",FNRD:"Fnord Software",FORA:"Fora",FORE:"Forefront Technology",FP:"Fujitsu",FPA:"WayTech Development",FUJI:"Fujitsu",FX:"Fuji Xerox",GCC:"GCC Technologies",GGSL:"Global Graphics Software",GMB:"Gretagmacbeth",GMG:"GMG",GOLD:"GoldStar Technology",GOOG:"Google",GPRT:"Giantprint",GTMB:"Gretagmacbeth",GVC:"WayTech Development",GW2K:"Sony",HCI:"HCI",HDM:"Heidelberger Druckmaschinen",HERM:"Hermes",HITA:"Hitachi America",HP:"Hewlett-Packard",HTC:"Hitachi",HiTi:"HiTi Digital",IBM:"IBM",IDNT:"Scitex",IEC:"Hewlett-Packard",IIYA:"Iiyama North America",IKEG:"Ikegami Electronics",IMAG:"Image Systems",IMI:"Ingram Micro",INTC:"Intel",INTL:"N/A (INTL)",INTR:"Intra Electronics",IOCO:"Iocomm International Technology",IPS:"InfoPrint Solutions Company",IRIS:"Scitex",ISL:"Ichikawa Soft Laboratory",ITNL:"N/A (ITNL)",IVM:"IVM",IWAT:"Iwatsu Electric",Idnt:"Scitex",Inca:"Inca Digital Printers",Iris:"Scitex",JPEG:"Joint Photographic Experts Group",JSFT:"Jetsoft Development",JVC:"JVC Information Products",KART:"Scitex",KFC:"KFC Computek Components",KLH:"KLH Computers",KMHD:"Konica Minolta",KNCA:"Konica",KODA:"Kodak",KYOC:"Kyocera",Kart:"Scitex",LCAG:"Leica",LCCD:"Leeds Colour",LDAK:"Left Dakota",LEAD:"Leading Technology",LEXM:"Lexmark International",LINK:"Link Computer",LINO:"Linotronic",LITE:"Lite-On",Leaf:"Leaf",Lino:"Linotronic",MAGC:"Mag Computronic",MAGI:"MAG Innovision",MANN:"Mannesmann",MICN:"Micron Technology",MICR:"Microtek",MICV:"Microvitec",MINO:"Minolta",MITS:"Mitsubishi Electronics America",MITs:"Mitsuba",MNLT:"Minolta",MODG:"Modgraph",MONI:"Monitronix",MONS:"Monaco Systems",MORS:"Morse Technology",MOTI:"Motive Systems",MSFT:"Microsoft",MUTO:"MUTOH INDUSTRIES",Mits:"Mitsubishi Electric",NANA:"NANAO",NEC:"NEC",NEXP:"NexPress Solutions",NISS:"Nissei Sangyo America",NKON:"Nikon",NONE:"none",OCE:"Oce Technologies",OCEC:"OceColor",OKI:"Oki",OKID:"Okidata",OKIP:"Okidata",OLIV:"Olivetti",OLYM:"Olympus",ONYX:"Onyx Graphics",OPTI:"Optiquest",PACK:"Packard Bell",PANA:"Matsushita Electric Industrial",PANT:"Pantone",PBN:"Packard Bell",PFU:"PFU",PHIL:"Philips Consumer Electronics",PNTX:"HOYA",POne:"Phase One A/S",PREM:"Premier Computer Innovations",PRIN:"Princeton Graphic Systems",PRIP:"Princeton Publishing Labs",QLUX:"Hong Kong",QMS:"QMS",QPCD:"QPcard AB",QUAD:"QuadLaser",QUME:"Qume",RADI:"Radius",RDDx:"Integrated Color Solutions",RDG:"Roland DG",REDM:"REDMS Group",RELI:"Relisys",RGMS:"Rolf Gierling Multitools",RICO:"Ricoh",RNLD:"Edmund Ronald",ROYA:"Royal",RPC:"Ricoh Printing Systems",RTL:"Royal Information Electronics",SAMP:"Sampo",SAMS:"Samsung",SANT:"Jaime Santana Pomares",SCIT:"Scitex",SCRN:"Dainippon Screen",SDP:"Scitex",SEC:"Samsung",SEIK:"Seiko Instruments",SEIk:"Seikosha",SGUY:"ScanGuy.com",SHAR:"Sharp Laboratories",SICC:"International Color Consortium",SONY:"Sony",SPCL:"SpectraCal",STAR:"Star",STC:"Sampo Technology",Scit:"Scitex",Sdp:"Scitex",Sony:"Sony",TALO:"Talon Technology",TAND:"Tandy",TATU:"Tatung",TAXA:"TAXAN America",TDS:"Tokyo Denshi Sekei",TECO:"TECO Information Systems",TEGR:"Tegra",TEKT:"Tektronix",TI:"Texas Instruments",TMKR:"TypeMaker",TOSB:"Toshiba",TOSH:"Toshiba",TOTK:"TOTOKU ELECTRIC",TRIU:"Triumph",TSBT:"Toshiba",TTX:"TTX Computer Products",TVM:"TVM Professional Monitor",TW:"TW Casper",ULSX:"Ulead Systems",UNIS:"Unisys",UTZF:"Utz Fehlau & Sohn",VARI:"Varityper",VIEW:"Viewsonic",VISL:"Visual communication",VIVO:"Vivo Mobile Communication",WANG:"Wang",WLBR:"Wilbur Imaging",WTG2:"Ware To Go",WYSE:"WYSE Technology",XERX:"Xerox",XRIT:"X-Rite",ZRAN:"Zoran",Zebr:"Zebra Technologies",appl:"Apple Computer",bICC:"basICColor",berg:"bergdesign",ceyd:"Integrated Color Solutions",clsp:"MacDermid ColorSpan",ds:"Dainippon Screen",dupn:"DuPont",ffei:"FujiFilm Electronic Imaging",flux:"FluxData",iris:"Scitex",kart:"Scitex",lcms:"Little CMS",lino:"Linotronic",none:"none",ob4d:"Erdt Systems",obic:"Medigraph",quby:"Qubyx Sarl",scit:"Scitex",scrn:"Dainippon Screen",sdp:"Scitex",siwi:"SIWI GRAFIKA",yxym:"YxyMaster"},Vr={scnr:"Scanner",mntr:"Monitor",prtr:"Printer",link:"Device Link",abst:"Abstract",spac:"Color Space Conversion Profile",nmcl:"Named Color",cenc:"ColorEncodingSpace profile",mid:"MultiplexIdentification profile",mlnk:"MultiplexLink profile",mvis:"MultiplexVisualization profile",nkpf:"Nikon Input Device Profile (NON-STANDARD!)"};Pt(le,"icc",[[4,Wn],[12,Vr],[40,Object.assign({},Wn,Vr)],[48,Wn],[80,Wn],[64,{0:"Perceptual",1:"Relative Colorimetric",2:"Saturation",3:"Absolute Colorimetric"}],["tech",{amd:"Active Matrix Display",crt:"Cathode Ray Tube Display",kpcd:"Photo CD",pmd:"Passive Matrix Display",dcam:"Digital Camera",dcpj:"Digital Cinema Projector",dmpc:"Digital Motion Picture Camera",dsub:"Dye Sublimation Printer",epho:"Electrophotographic Printer",esta:"Electrostatic Printer",flex:"Flexography",fprn:"Film Writer",fscn:"Film Scanner",grav:"Gravure",ijet:"Ink Jet Printer",imgs:"Photo Image Setter",mpfr:"Motion Picture Film Recorder",mpfs:"Motion Picture Film Scanner",offs:"Offset Lithography",pjtv:"Projection Television",rpho:"Photographic Paper Printer",rscn:"Reflective Scanner",silk:"Silkscreen",twax:"Thermal Wax Printer",vidc:"Video Camera",vidm:"Video Monitor"}]]);class jn extends oe{static canHandle(t,e,i){return t.getUint8(e+1)===237&&t.getString(e+4,9)==="Photoshop"&&this.containsIptc8bim(t,e,i)!==void 0}static headerLength(t,e,i){let r,s=this.containsIptc8bim(t,e,i);if(s!==void 0)return r=t.getUint8(e+s+7),r%2!=0&&(r+=1),r===0&&(r=4),s+8+r}static containsIptc8bim(t,e,i){for(let r=0;r<i;r++)if(this.isIptcSegmentHead(t,e+r))return r}static isIptcSegmentHead(t,e){return t.getUint8(e)===56&&t.getUint32(e)===943868237&&t.getUint16(e+4)===1028}parse(){let{raw:t}=this,e=this.chunk.byteLength-1,i=!1;for(let r=0;r<e;r++)if(this.chunk.getUint8(r)===28&&this.chunk.getUint8(r+1)===2){i=!0;let s=this.chunk.getUint16(r+3),o=this.chunk.getUint8(r+2),a=this.chunk.getLatin1String(r+5,s);t.set(o,this.pluralizeValue(t.get(o),a)),r+=4+s}else if(i)break;return this.translate(),this.output}pluralizeValue(t,e){return t!==void 0?t instanceof Array?(t.push(e),t):[t,e]:e}}st(jn,"type","iptc"),st(jn,"translateValues",!1),st(jn,"reviveValues",!1),Ct.set("iptc",jn),Pt(Rt,"iptc",[[0,"ApplicationRecordVersion"],[3,"ObjectTypeReference"],[4,"ObjectAttributeReference"],[5,"ObjectName"],[7,"EditStatus"],[8,"EditorialUpdate"],[10,"Urgency"],[12,"SubjectReference"],[15,"Category"],[20,"SupplementalCategories"],[22,"FixtureIdentifier"],[25,"Keywords"],[26,"ContentLocationCode"],[27,"ContentLocationName"],[30,"ReleaseDate"],[35,"ReleaseTime"],[37,"ExpirationDate"],[38,"ExpirationTime"],[40,"SpecialInstructions"],[42,"ActionAdvised"],[45,"ReferenceService"],[47,"ReferenceDate"],[50,"ReferenceNumber"],[55,"DateCreated"],[60,"TimeCreated"],[62,"DigitalCreationDate"],[63,"DigitalCreationTime"],[65,"OriginatingProgram"],[70,"ProgramVersion"],[75,"ObjectCycle"],[80,"Byline"],[85,"BylineTitle"],[90,"City"],[92,"Sublocation"],[95,"State"],[100,"CountryCode"],[101,"Country"],[103,"OriginalTransmissionReference"],[105,"Headline"],[110,"Credit"],[115,"Source"],[116,"CopyrightNotice"],[118,"Contact"],[120,"Caption"],[121,"LocalCaption"],[122,"Writer"],[125,"RasterizedCaption"],[130,"ImageType"],[131,"ImageOrientation"],[135,"LanguageIdentifier"],[150,"AudioType"],[151,"AudioSamplingRate"],[152,"AudioSamplingResolution"],[153,"AudioDuration"],[154,"AudioOutcue"],[184,"JobID"],[185,"MasterDocumentID"],[186,"ShortDocumentID"],[187,"UniqueDocumentID"],[188,"OwnerID"],[200,"ObjectPreviewFileFormat"],[201,"ObjectPreviewFileVersion"],[202,"ObjectPreviewData"],[221,"Prefs"],[225,"ClassifyState"],[228,"SimilarityIndex"],[230,"DocumentNotes"],[231,"DocumentHistory"],[232,"ExifCameraInfo"],[255,"CatalogSets"]]),Pt(le,"iptc",[[10,{0:"0 (reserved)",1:"1 (most urgent)",2:"2",3:"3",4:"4",5:"5 (normal urgency)",6:"6",7:"7",8:"8 (least urgent)",9:"9 (user-defined priority)"}],[75,{a:"Morning",b:"Both Morning and Evening",p:"Evening"}],[131,{L:"Landscape",P:"Portrait",S:"Square"}]]);function Oa(n,t,e){const i=e.type===1||e.type===2||e.type===7?1:e.type===3?2:e.type===4||e.type===9||e.type===11||e.type===13?4:e.type===5||e.type===10||e.type===12?8:0;if(i===0)throw new Error(`Unsupported TIFF field type: ${e.type}`);const s=e.count*i<=4?e.valueFieldOffset:n.getUint32(e.valueFieldOffset,t),o=[];for(let a=0;a<e.count;a++){const l=s+a*i;if(l<0||l+i>n.byteLength)throw new Error("Invalid TIFF field offset");if(e.type===1)o.push(n.getUint8(l));else if(e.type===2||e.type===7)o.push(n.getUint8(l));else if(e.type===3)o.push(n.getUint16(l,t));else if(e.type===5){const u=n.getUint32(l+4,t);o.push(u?n.getUint32(l,t)/u:0)}else if(e.type===13)o.push(n.getUint32(l,t));else if(e.type===9)o.push(n.getInt32(l,t));else if(e.type===10){const u=n.getInt32(l+4,t);o.push(u?n.getInt32(l,t)/u:0)}else e.type===11?o.push(n.getFloat32(l,t)):e.type===12?o.push(n.getFloat64(l,t)):o.push(n.getUint32(l,t))}return o}function Xr(n,t,e){if(!Number.isFinite(e)||e<=0||e+2>n.byteLength)throw new Error("Invalid TIFF IFD offset");const i=n.getUint16(e,t);if(e+2+i*12+4>n.byteLength)throw new Error("Corrupt TIFF IFD");const s=new Map;for(let o=0;o<i;o++){const a=e+2+o*12,l=n.getUint16(a,t),u=n.getUint16(a+2,t),h=n.getUint32(a+4,t);s.set(l,{type:u,count:h,valueFieldOffset:a+8})}return s}function Wt(n,t,e,i){const r=e.get(i);return r?Oa(n,t,r):[]}const ys=(...n)=>{for(const t of n){if(Array.isArray(t)||ArrayBuffer.isView(t)){const i=Array.from(t).map(Number).filter(r=>Number.isFinite(r)&&r>0);if(i.length)return Math.max(...i);continue}const e=Number(t);if(Number.isFinite(e)&&e>0)return e}return null},Va=n=>{const e=(Array.isArray(n)||ArrayBuffer.isView(n)?Array.from(n).map(Number):[Number(n)]).map(i=>Number.isFinite(i)?i:0);return e.length>=3?[e[0]||0,e[1]||0,e[1]||0,e[2]||0]:e.length===1?[e[0]||0,e[0]||0,e[0]||0,e[0]||0]:[0,0,0,0]};function Xa(n,t={}){var A,F,C,M,S,I;if(n.byteLength<8)return null;const e=new DataView(n.buffer,n.byteOffset,n.byteLength),i=e.getUint16(0,!1),r=i===18761;if(!r&&i!==19789||e.getUint16(2,r)!==42)return null;const s=e.getUint32(4,r);let o;try{o=Xr(e,r,s)}catch{return null}const a=new Set([s]);for(const N of Wt(e,r,o,330))N>0&&a.add(N);let l=null;for(const N of a)try{const E=N===s?o:Xr(e,r,N),X=Wt(e,r,E,256)[0]||0,U=Wt(e,r,E,257)[0]||0,L=Wt(e,r,E,258),B=Wt(e,r,E,277)[0]||L.length||1,O=L.length===1?new Array(B).fill(L[0]):L,z=Wt(e,r,E,259)[0]||1,J=Wt(e,r,E,262)[0]||0,V=Wt(e,r,E,284)[0]||1,q=Wt(e,r,E,339),W=q.length===0?new Array(B).fill(1):q.length===1?new Array(B).fill(q[0]):q,tt=Wt(e,r,E,273),at=Wt(e,r,E,279),D=Wt(e,r,E,278)[0]||U,Z=W.slice(0,B).every(T=>T===0||T===1);if(!(J===34892&&z===1&&V===1&&X>0&&U>0&&B>=3&&O.length>=B&&O.slice(0,B).every(T=>T===16)&&Z&&tt.length>0&&tt.length===at.length))continue;(!l||X*U>l.width*l.height)&&(l={entries:E,width:X,height:U,samplesPerPixel:B,rowsPerStrip:D,stripOffsets:tt,stripByteCounts:at,bitsPerSample:O,photometric:J})}catch{continue}if(!l)return null;const{entries:u,width:h,height:c,samplesPerPixel:g,rowsPerStrip:d,stripOffsets:f,stripByteCounts:p,bitsPerSample:m,photometric:y}=l,x=h*c;if(!Number.isSafeInteger(x)||x<=0)return null;const b=new Uint16Array(x*3);for(let N=0;N<f.length;N++){const E=N*d;if(E>=c)break;const X=Math.min(d,c-E),U=X*h*g*2,L=f[N],B=p[N];if(L<0||B<U||L+U>n.byteLength)throw new Error("Invalid LinearRaw DNG strip bounds");const O=E*h*3;if(g===3&&r&&!(n.byteOffset+L&1)){b.set(new Uint16Array(n.buffer,n.byteOffset+L,X*h*3),O);continue}const z=new DataView(n.buffer,n.byteOffset+L,U);let J=0,V=O;for(let q=0;q<X*h;q++)b[V++]=z.getUint16(J,r),b[V++]=z.getUint16(J+2,r),b[V++]=z.getUint16(J+4,r),J+=g*2}const _=Wt(e,r,u,50714),k=((F=(A=t==null?void 0:t.color_data)==null?void 0:A.dng_levels)==null?void 0:F.dng_cblack)||((C=t==null?void 0:t.color_data)==null?void 0:C.cblack_rawpy_style)||(t==null?void 0:t.black_level_per_channel)||(t==null?void 0:t.cblack),v=Va(_.length?_:k),w=Wt(e,r,u,50717),P=ys(w,(S=(M=t==null?void 0:t.color_data)==null?void 0:M.dng_levels)==null?void 0:S.dng_whitelevel,(I=t==null?void 0:t.color_data)==null?void 0:I.maximum,t==null?void 0:t.white_level)||65535;return{data:b,width:h,height:c,bayerPattern:"",blackLevels:v,whiteLevel:P,metadata:{...t,format:"DNG_LINEAR_RAW_RGB",description:"Uncompressed DNG LinearRaw RGB",linearRawDngDecoder:!0,bitsPerSample:m.slice(0,g),samplesPerPixel:g,photometric:y},isThreePlane:!0,threePlaneTransfer:"linear"}}async function Ga(n){const t=Xa(n);if(!t)return null;let e={};const i=await Qi(),r=new i;try{await r.open(n,{}),e=await r.metadata(!0)}catch(s){console.warn("LinearRaw DNG metadata enrichment failed",s)}finally{try{r.delete?r.delete():r.close()}catch{}}try{const s=await gs.parse(n.buffer);s&&(e={...e,...s})}catch(s){console.warn("exifr parsing failed for LinearRaw DNG",s)}return t.metadata={...e,...t.metadata},t}const za=async n=>{var r,s,o,a,l,u,h,c,g,d,f;const t=await Ga(n);if(t)return t;const e=await Qi(),i=new e;try{if(await i.open(n,{}),typeof i.getRawImage!="function")throw new Error("WASM mismatch");let p={};try{p=await i.metadata(!0)}catch(S){console.warn("Metadata error before raw extraction",S)}const m=i.getRawImage(),y=m.data instanceof Uint16Array?m.data:new Uint16Array(m.data);let x={...p};try{const S=await gs.parse(n.buffer);S&&(x={...x,...S})}catch(S){console.warn("exifr parsing failed for RAW buffer",S)}const b=p.filters||((r=p.idata)==null?void 0:r.filters)||0,_=p.colors||((s=p.idata)==null?void 0:s.colors)||0,k=b===0&&_===3,v=b===9;let w=[0,0,0,0],P=null,A;if(i.getBlackLevels)try{const S=i.getBlackLevels();P=S,A=wr(S)||wr((o=p.color_data)==null?void 0:o.black_level_model)||void 0;const I=$e(A==null?void 0:A.siteBaseLevels);I&&(w=I)}catch(S){console.warn("getBlackLevels binding failed",S)}if(!A){const S=Number((P==null?void 0:P.black)??((a=p.color_data)==null?void 0:a.black)??0)||0,I=$e((P==null?void 0:P.cblack)||((l=p.color_data)==null?void 0:l.cblack_rawpy_style)||p.black_level_per_channel||p.cblack||((u=p.color)==null?void 0:u.cblack));if(I){const N=String(m.bayerPattern||p.cfa_pattern||"RGGB").toUpperCase();let E=0;w=[0,1,2,3].map(X=>{const U=N[X],L=U==="R"?0:U==="B"?2:++E===1?1:3;return Math.max(0,S+I[L])})}else w=[S,S,S,S]}const F=Number.isFinite(Number(m.bits))&&Number(m.bits)>0?Math.pow(2,Number(m.bits))-1:null,C=ys(p.white_level,(c=(h=p.color_data)==null?void 0:h.dng_levels)==null?void 0:c.dng_whitelevel,(g=p.color_data)==null?void 0:g.maximum,P==null?void 0:P.maximum,(d=p.color_data)==null?void 0:d.fmaximum,(f=p.color_data)==null?void 0:f.data_maximum,F)||16383,M={data:y,width:m.width,height:m.height,bayerPattern:m.bayerPattern||"",blackLevels:w,blackLevelModel:A,whiteLevel:C,metadata:x,isThreePlane:k,threePlaneTransfer:k?"linear":void 0,isXTrans:v};return fa(M,da(M,x)),M}finally{i.delete?i.delete():i.close()}};async function Ya(n){if(ia(n)){const i=await la(n);if(!i)throw new Error("Sony cRAW HQ decoder did not return image data.");return i.rawImageData}const e=new Uint8Array(n);return za(e)}function Wa(n,t,e){const i=xs(n,e),r=Math.floor(t.x),s=Math.floor(t.y),o=Math.floor(t.w),a=Math.floor(t.h),l=new Uint16Array(o*a);for(let u=0;u<a;u++){const h=s+u,c=u*o;for(let g=0;g<o;g++)l[c+g]=i(r+g,h)}return{data:l,width:o,height:a}}function ja(n,t,e){const i=Ha(n,t,e);if(!i)return null;const r=n.width,s=n.height,o=new Uint16Array(r*s);for(let a=0;a<s;a++){const l=a*r;for(let u=0;u<r;u++)o[l+u]=i(u,a)}return{kind:"u16-mono",data:o,width:r,height:s}}function Ha(n,t,e){return t.renderMode==="advanced-zero-dep"&&t.advancedZeroDep?xs(n,t,e):t.renderMode==="zero-dependency"?Qa(n,t,e):null}function xs(n,t,e){if(!t.advancedZeroDep)throw new Error("Unmixing settings not found in DisplaySettings.");const{bg:i,fg:r}=t.advancedZeroDep,s=bs(e,t.advancedZeroDep.bl),{data:o,width:a,whiteLevel:l}=n,u=i.map((f,p)=>Math.max(0,f-s[p])),h=r.map((f,p)=>Math.max(0,f-s[p])),c=(u[1]+u[3])/2,g=(h[1]+h[3])/2,d=Math.pow(2,t.exposure);return(f,p)=>{if(f<0||p<0||f>=a||p>=n.height)return 0;const m=p%2,y=f%2;let x=0;!m&&!y?x=0:!m&&y?x=1:m&&!y?x=3:x=2;const b=o[p*a+f],_=s[x],k=u[x],v=h[x],w=Math.max(b-_,0),P=v-k||1e-9,A=(w-k)/P;let F;return A<0?F=w*(c/Math.max(k,1e-9)):A>1?F=w*(g/Math.max(v,1e-9)):F=(1-A)*c+A*g,F*=d,Math.max(0,Math.min(65535,Math.round(F)))}}function Qa(n,t,e){const{data:i,width:r,height:s,whiteLevel:o}=n,a=bs(e,t.blackLevel||[0,0,0,0]),l=Ka(n.bayerPattern),u=t.wbGains?t.wbGains[0]:1,h=t.wbGains?t.wbGains[1]:1,c=Math.pow(2,t.exposure||0);return(g,d)=>{if(g<0||d<0||g>=r||d>=s)return 0;const f=$a(g,d),p=Ja(l,g,d),m=i[d*r+g],y=a[f];let x=(m-y)/Math.max(1,o-y);return x=Math.max(0,Math.min(1,x)),x*=c,p==="R"?x*=u:p==="B"&&(x*=h),qa(x)}}function bs(n,t){if(typeof n=="number"&&Number.isFinite(n)){const e=Math.max(0,n);return[e,e,e,e]}return Array.isArray(n)&&n.length===4?[Math.max(0,n[0]??0),Math.max(0,n[1]??0),Math.max(0,n[2]??0),Math.max(0,n[3]??0)]:[Math.max(0,t[0]??0),Math.max(0,t[1]??0),Math.max(0,t[2]??0),Math.max(0,t[3]??0)]}function qa(n){return Math.max(0,Math.min(65535,Math.round(Math.max(0,Math.min(1,n))*65535)))}function Ka(n){const t=Je(n);if(!t)throw new Error("Cannot process RAW mosaic without a valid Bayer CFA pattern.");return t}function $a(n,t){return ns(n,t)}function Ja(n,t,e){return pa(n,t,e)}const ft=(n,t=0)=>({real:n,imag:t}),gn=(n,t)=>({real:n.real+t.real,imag:n.imag+t.imag}),pn=(n,t)=>({real:n.real-t.real,imag:n.imag-t.imag}),It=(n,t)=>({real:n.real*t.real-n.imag*t.imag,imag:n.real*t.imag+n.imag*t.real}),Ue=(n,t)=>{const e=t.real*t.real+t.imag*t.imag;return e===0?ft(0):{real:(n.real*t.real+n.imag*t.imag)/e,imag:(n.imag*t.real-n.real*t.imag)/e}},Ke=n=>Math.hypot(n.real,n.imag),_s=n=>{const t=Ke(n);if(t===0)return ft(0);const e=Math.sqrt(t),i=Math.atan2(n.imag,n.real);return ft(e*Math.cos(i/2),e*Math.sin(i/2))};function Za(n,t){const e=n.length-1;if(e<0)return{p:ft(0),dp:ft(0),d2p:ft(0)};let i=ft(n[e].real,n[e].imag),r=ft(0),s=ft(0);for(let o=e-1;o>=0;o--)s=gn(It(r,ft(2)),It(t,s)),r=gn(i,It(t,r)),i=gn(ft(n[o].real,n[o].imag),It(t,i));return{p:i,dp:r,d2p:s}}function Wi(n,t,e=80){const r=n.length-1;if(r<=0)return{root:t,iterations:0};if(r===1)return{root:Ue(It(n[0],ft(-1)),n[1]),iterations:0};let s=ft(t.real,t.imag);for(let o=0;o<e;o++){const{p:a,dp:l,d2p:u}=Za(n,s);if(Ke(a)<1e-14)return{root:s,iterations:o};const h=Ue(l,a),c=It(h,h),g=pn(c,Ue(u,a)),d=ft(r),f=ft(r-1),p=pn(It(d,g),It(h,h)),m=_s(It(f,p)),y=gn(h,m),x=pn(h,m),b=Ke(y)>Ke(x)?y:x;if(Ke(b)<1e-14)return{root:s,iterations:o};const _=Ue(d,b),k=pn(s,_);if(Ke(_)<1e-14*Ke(k))return{root:k,iterations:o+1};s=k}return{root:s,iterations:e}}function to(n,t){const e=n.length-1;if(e<=0)return[ft(0)];if(e===1)return[ft(n[0].real,n[0].imag)];const i=new Array(e);i[e-1]=ft(n[e].real,n[e].imag);for(let r=e-2;r>=0;r--){const s=ft(n[r+1].real,n[r+1].imag),o=i[r+1];i[r]=gn(s,It(t,o))}return i}function eo(n){const t=n.length-1;if(t<=0)return[];if(t===1)return[Ue(It(n[0],ft(-1)),n[1])];const e=[];let i=n.map(s=>ft(s.real,s.imag)),r=t*5;for(;i.length>2&&r-- >0;){const s=ft(.3+Math.random()*.7,.3+Math.random()*.7),{root:o}=Wi(i,s,100),a=Wi(n,o,20);e.push(a.root);const l=to(i,o);if(l.length>=i.length){console.warn("polyDeflate did not reduce degree, breaking");break}i=l}if(i.length===2)e.push(Ue(It(i[0],ft(-1)),i[1]));else if(i.length===3){const s=i[2],o=i[1],a=i[0],l=pn(It(o,o),It(It(ft(4),s),a)),u=_s(l),h=It(ft(2),s),c=Ue(pn(It(ft(-1),o),u),h),g=Ue(gn(It(ft(-1),o),u),h);e.push(c,g)}return e}function no(n,t,e){const i=[ft(n),ft(-1),ft(n*t),ft(0),ft(n*e)],r=eo(i);if(r.length===0)return console.warn("laguerreSmallestPositiveRoot: no roots found"),n;let s=1/0,o=!1;for(const l of r)Math.abs(l.imag)<1e-10&&l.real>0&&l.real<s&&(s=l.real,o=!0);return o?Wi(i,ft(s,0),20).root.real:(console.warn("laguerreSmallestPositiveRoot: no positive real root found"),n)}function io(n,t,e){if(Math.abs(t)<1e-10&&Math.abs(e)<1e-10)return n;if(n<1e-10)return 0;if(Math.abs(e)<1e-10){const i=-1/(t*n),r=1/t,s=i*i-4*r;if(s<0)return n;const o=Math.sqrt(s),a=-.5*(i+Math.sign(i)*o),l=a,u=r/a;return l>0&&u>0?Math.min(l,u):l>0?l:u>0?u:n}try{return no(n,t,e)}catch(i){return console.error("laguerreSmallestPositiveRoot failed:",i),n}}function ro(n,t,e,i,r){const s=n.x-t.x,o=n.y-t.y,a=Math.hypot(s,o)/Math.max(1e-12,e);if(a<1e-12)return{x:n.x,y:n.y};const l=a*a,u=1+(i+r*l)*l;return{x:s/u+t.x,y:o/u+t.y}}function so(n,t,e,i,r){const s=n.x-t.x,o=n.y-t.y,a=Math.hypot(s,o)/Math.max(1e-12,e);if(a<1e-12)return{x:t.x,y:t.y};const u=io(a,i,r)/a;return{x:t.x+s*u,y:t.y+o*u}}function Ut(n,t){const e=ro(n,{x:t.principalX,y:t.principalY},t.radiusNorm,t.k1,t.k2);return{x:e.x+(t.correctedOffsetX??0),y:e.y+(t.correctedOffsetY??0)}}function Fe(n,t){const e={x:n.x-(t.correctedOffsetX??0),y:n.y-(t.correctedOffsetY??0)};return so(e,{x:t.principalX,y:t.principalY},t.radiusNorm,t.k1,t.k2)}const ao=`
attribute vec2 a_position;

void main() {
  gl_Position = vec4(a_position, 0.0, 1.0);
}
`,oo=`
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
`,lo=`
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
`;class co{constructor(){bt(this,"canvas",null);bt(this,"gl",null);bt(this,"blurProgram",null);bt(this,"sobelProgram",null);bt(this,"positionBuffer",null);bt(this,"blurUniforms",null);bt(this,"sobelUniforms",null);bt(this,"resources",null);bt(this,"initialized",!1);bt(this,"unavailable",!1);bt(this,"maxTextureSize",0)}compute(t,e,i){if(!this.initialized&&!this.init())return null;const r=this.gl,s=this.blurProgram,o=this.sobelProgram,a=this.blurUniforms,l=this.sobelUniforms;if(!r||!s||!o||!a||!l||!this.positionBuffer||!this.canvas||e<=2||i<=2||e>this.maxTextureSize||i>this.maxTextureSize)return null;const u=this.ensureResources(e,i);if(!u)return null;this.canvas.width=e,this.canvas.height=i,r.viewport(0,0,e,i),r.disable(r.BLEND),r.pixelStorei(r.UNPACK_ALIGNMENT,1),r.pixelStorei(r.PACK_ALIGNMENT,1),r.bindBuffer(r.ARRAY_BUFFER,this.positionBuffer),r.activeTexture(r.TEXTURE0),r.bindTexture(r.TEXTURE_2D,u.sourceTexture),r.texImage2D(r.TEXTURE_2D,0,r.LUMINANCE,e,i,0,r.LUMINANCE,r.UNSIGNED_BYTE,t),r.useProgram(s),r.enableVertexAttribArray(0),r.vertexAttribPointer(0,2,r.FLOAT,!1,0,0),r.uniform2f(a.size,e,i),r.uniform1i(a.source,0),r.bindFramebuffer(r.FRAMEBUFFER,u.blurFramebuffer),r.bindTexture(r.TEXTURE_2D,u.sourceTexture),r.drawArrays(r.TRIANGLES,0,6);const h=new Uint8Array(e*i*4);r.readPixels(0,0,e,i,r.RGBA,r.UNSIGNED_BYTE,h),r.useProgram(o),r.uniform2f(l.size,e,i),r.uniform1i(l.blurred,0),r.bindFramebuffer(r.FRAMEBUFFER,u.sobelFramebuffer),r.bindTexture(r.TEXTURE_2D,u.blurTexture),r.drawArrays(r.TRIANGLES,0,6);const c=new Uint8Array(e*i*4);r.readPixels(0,0,e,i,r.RGBA,r.UNSIGNED_BYTE,c),r.disableVertexAttribArray(0),r.bindFramebuffer(r.FRAMEBUFFER,null),r.bindBuffer(r.ARRAY_BUFFER,null),r.bindTexture(r.TEXTURE_2D,null);const g=new Uint8Array(e*i);for(let m=0,y=0;m<g.length;m++,y+=4)g[m]=h[y];const d=new Float32Array(e*i),f=new Float32Array(e*i),p=new Float32Array(e*i);for(let m=0,y=0;m<d.length;m++,y+=4){const x=(c[y]|c[y+1]<<8)-32768,b=(c[y+2]|c[y+3]<<8)-32768;d[m]=x,f[m]=b,p[m]=Math.sqrt(x*x+b*b)}return{blurredGray:g,gx:d,gy:f,magnitude:p}}init(){if(this.initialized&&this.gl&&this.blurProgram&&this.sobelProgram)return!0;if(this.unavailable)return!1;const t=this.createCanvas();if(!t)return this.unavailable=!0,!1;const e=t.getContext("webgl",{alpha:!1,antialias:!1,depth:!1,stencil:!1,premultipliedAlpha:!1,preserveDrawingBuffer:!1});if(!e)return this.unavailable=!0,!1;const i=this.compileShader(e,e.VERTEX_SHADER,ao),r=this.compileShader(e,e.FRAGMENT_SHADER,oo),s=this.compileShader(e,e.FRAGMENT_SHADER,lo);if(!i||!r||!s)return i&&e.deleteShader(i),r&&e.deleteShader(r),s&&e.deleteShader(s),this.unavailable=!0,!1;const o=this.createProgram(e,i,r),a=this.createProgram(e,i,s);if(e.deleteShader(i),e.deleteShader(r),e.deleteShader(s),!o||!a)return o&&e.deleteProgram(o),a&&e.deleteProgram(a),this.unavailable=!0,!1;const l=e.createBuffer();return l?(e.bindBuffer(e.ARRAY_BUFFER,l),e.bufferData(e.ARRAY_BUFFER,new Float32Array([-1,-1,1,-1,-1,1,-1,1,1,-1,1,1]),e.STATIC_DRAW),e.bindBuffer(e.ARRAY_BUFFER,null),this.canvas=t,this.gl=e,this.blurProgram=o,this.sobelProgram=a,this.positionBuffer=l,this.blurUniforms={source:e.getUniformLocation(o,"u_source"),size:e.getUniformLocation(o,"u_size")},this.sobelUniforms={blurred:e.getUniformLocation(a,"u_blurred"),size:e.getUniformLocation(a,"u_size")},this.maxTextureSize=Number(e.getParameter(e.MAX_TEXTURE_SIZE)||0),this.initialized=!0,!0):(e.deleteProgram(o),e.deleteProgram(a),this.unavailable=!0,!1)}createCanvas(){return typeof OffscreenCanvas<"u"?new OffscreenCanvas(1,1):typeof document<"u"?document.createElement("canvas"):null}ensureResources(t,e){const i=this.gl;if(!i)return null;if(this.resources&&this.resources.width===t&&this.resources.height===e)return this.resources;this.disposeResources();const r=this.createTexture(i.LUMINANCE,t,e,i.LUMINANCE,i.UNSIGNED_BYTE,null),s=this.createTexture(i.RGBA,t,e,i.RGBA,i.UNSIGNED_BYTE,null),o=this.createTexture(i.RGBA,t,e,i.RGBA,i.UNSIGNED_BYTE,null),a=this.createFramebuffer(s),l=this.createFramebuffer(o);return!r||!s||!o||!a||!l?(r&&i.deleteTexture(r),s&&i.deleteTexture(s),o&&i.deleteTexture(o),a&&i.deleteFramebuffer(a),l&&i.deleteFramebuffer(l),null):(this.resources={width:t,height:e,sourceTexture:r,blurTexture:s,sobelTexture:o,blurFramebuffer:a,sobelFramebuffer:l},this.resources)}createTexture(t,e,i,r,s,o){const a=this.gl;if(!a)return null;const l=a.createTexture();return l?(a.bindTexture(a.TEXTURE_2D,l),a.texParameteri(a.TEXTURE_2D,a.TEXTURE_MIN_FILTER,a.NEAREST),a.texParameteri(a.TEXTURE_2D,a.TEXTURE_MAG_FILTER,a.NEAREST),a.texParameteri(a.TEXTURE_2D,a.TEXTURE_WRAP_S,a.CLAMP_TO_EDGE),a.texParameteri(a.TEXTURE_2D,a.TEXTURE_WRAP_T,a.CLAMP_TO_EDGE),a.texImage2D(a.TEXTURE_2D,0,t,e,i,0,r,s,o),a.bindTexture(a.TEXTURE_2D,null),l):null}createFramebuffer(t){const e=this.gl;if(!e||!t)return null;const i=e.createFramebuffer();if(!i)return null;e.bindFramebuffer(e.FRAMEBUFFER,i),e.framebufferTexture2D(e.FRAMEBUFFER,e.COLOR_ATTACHMENT0,e.TEXTURE_2D,t,0);const r=e.checkFramebufferStatus(e.FRAMEBUFFER);return e.bindFramebuffer(e.FRAMEBUFFER,null),r!==e.FRAMEBUFFER_COMPLETE?(e.deleteFramebuffer(i),null):i}compileShader(t,e,i){const r=t.createShader(e);return r?(t.shaderSource(r,i),t.compileShader(r),t.getShaderParameter(r,t.COMPILE_STATUS)?r:(console.error("[SFR Auto Detect WebGL] shader compile failed",t.getShaderInfoLog(r)),t.deleteShader(r),null)):null}createProgram(t,e,i){const r=t.createProgram();return r?(t.attachShader(r,e),t.attachShader(r,i),t.bindAttribLocation(r,0,"a_position"),t.linkProgram(r),t.getProgramParameter(r,t.LINK_STATUS)?r:(console.error("[SFR Auto Detect WebGL] program link failed",t.getProgramInfoLog(r)),t.deleteProgram(r),null)):null}disposeResources(){const t=this.gl,e=this.resources;if(!t||!e){this.resources=null;return}t.deleteTexture(e.sourceTexture),t.deleteTexture(e.blurTexture),t.deleteTexture(e.sobelTexture),t.deleteFramebuffer(e.blurFramebuffer),t.deleteFramebuffer(e.sobelFramebuffer),this.resources=null}}const uo=new co,dn=(n,t)=>{const e=Je(n);if(!e)throw new Error(`Bayer CFA pattern is unresolved for ${t}.`);return e},ji={gradientPercentiles:[.82,.88,.92,.95,.98,.995],downsampleMaxSide:1600,minComponentAreaRatio:15e-6,maxComponentAreaRatio:.35,minComponentAreaPx:20,minEdgePoints:24,extentQuantileLow:.02,extentQuantileHigh:.98,cornerTrimRatio:.18,minSpanPx:8,maxAspectRatio:2,bandScale:.16,bandMinPx:1.75,bandMaxPx:14,minPointContrast:6,minSidePoints:3,minCoverageRatio:.15,minCenterCoverageRatio:.2,filterBlockPurity:!0,innerPurityStdScale:1.5,outerMeanSpreadLimit:51,minAxisDot:.6,residualLimitFloor:.01,residualLimitScale:.25,minQuadArea:48,minSideLength:10,minOuterContrast:5,sampleHalfWidthRatio:.25};function ho(n,t,e,i,r,s){const o=n.width,a=n.height,l=dn(n.bayerPattern,"corrected RAW SFR sampling"),u=[],h=[],c=s!=null&&s.correctedRect?Lt*2:Lt,g=Math.max(1,Math.min(r,c)),d=e.p2.x-e.p1.x,f=e.p2.y-e.p1.y,p=Math.hypot(d,f);if(!Number.isFinite(p)||p<=1e-6)return null;const m=d/p,y=f/p,x=-y,b=m,_={x:(e.p1.x+e.p2.x)*.5,y:(e.p1.y+e.p2.y)*.5},k=s!=null&&s.correctedRect?Mt(s.correctedRect,o,a):Mt(Vt(_e(e,c*4+2)??[e.p1,e.p2],2),o,a);if(!k)return null;const v=(s==null?void 0:s.correctedScanlinesOverride)??(s!=null&&s.distortedRect?lr(Mt(s.distortedRect,o,a)??s.distortedRect,t,o,a):Us(k,e,Math.max(1,i),g*4+.5,o,a));if(!v||v.size===0)return null;const w=Bs(v,t,o,a);if(w.size===0)return null;const P=!cr(t);for(const[F,C]of w){if(F<0||F>=a)continue;const M=F*o;for(let S=C.start;S<=C.end;S++){if(S<0||S>=o||!vt(S,F,l,s==null?void 0:s.greenPhase))continue;const I={x:S,y:F},N=Ut(I,t);if(!Number.isFinite(N.x)||!Number.isFinite(N.y)||Math.round(N.x)<0||Math.round(N.x)>=o||Math.round(N.y)<0||Math.round(N.y)>=a)continue;const E=N.x-_.x,X=N.y-_.y,U=E*m+X*y;let L=E*x+X*b;if(P){const B=Ds(U,m,y,_,I,t);if(!B)continue;const O=.5*(B.a+B.b),z=Kt(B.a,m,y,_,t),J=Kt(O,m,y,_,t),V=Kt(B.b,m,y,_,t),q=hr({x:B.a,y:Math.hypot(z.x-I.x,z.y-I.y)},{x:O,y:Math.hypot(J.x-I.x,J.y-I.y)},{x:B.b,y:Math.hypot(V.x-I.x,V.y-I.y)});if(!Number.isFinite(q))continue;const W=ur(q,m,y,_,t),tt=Math.hypot(W.x,W.y);if(!Number.isFinite(tt)||tt<=1e-9)continue;const at=W.x/tt,Z=-(W.y/tt),H=at,T=Kt(q,m,y,_,t);L=(I.x-T.x)*Z+(I.y-T.y)*H}!Number.isFinite(U)||Math.abs(U)>Math.max(1,i)||!Number.isFinite(L)||Math.abs(L)>g||(u.push(L),h.push(Math.max(0,n.data[M+S]-Mn(s==null?void 0:s.blackLevel,S,F))))}}if(u.length<8)return null;const A=Math.abs(d)>=Math.abs(f)?1:2;return s!=null&&s.forceLegacyModel?Rn(u,h,A,c):wn(u,h,A,c)}function ir(n,t){const e=n.length;let i=0,r=0,s=0,o=0;for(let l=0;l<e;l++)i+=n[l],r+=t[l],s+=n[l]*t[l],o+=n[l]*n[l];const a=e*o-i*i;return a===0?{slope:0,intercept:0}:{slope:(e*s-i*r)/a,intercept:(r*o-i*s)/a}}function fo(n,t){const e=n.length,i=new Array(e).fill(0),r=2*t;for(let s=0;s<e;s++){const o=s>0?n[s-1]:n[0],a=s<e-1?n[s+1]:n[e-1];i[s]=(a-o)/r}return i}const te=-1e7,yn=13,ae=512,Et=8,fi=1/Et,Lt=28,po=[[0,0,0,0,0,-.085714285714286,.342857142857143,.485714285714286,.342857142857143,-.085714285714286,0,0,0,0,0],[0,0,0,0,-.095238095238095,.142857142857143,.285714285714286,.333333333333333,.285714285714286,.142857142857143,-.095238095238095,0,0,0,0],[0,0,0,-.090909090909091,.060606060606061,.168831168831169,.233766233766234,.255411255411255,.233766233766234,.168831168831169,.060606060606061,-.090909090909091,0,0,0],[0,0,-.083916083916084,.020979020979021,.102564102564103,.160839160839161,.195804195804196,.207459207459208,.195804195804196,.160839160839161,.102564102564103,.020979020979021,-.083916083916084,0,0],[0,-.076923076923077,0,.062937062937063,.111888111888112,.146853146853147,.167832167832168,.174825174825175,.167832167832168,.146853146853147,.111888111888112,.062937062937063,0,-.076923076923077,0],[-.070588235294118,-.011764705882353,.038009049773756,.078733031674208,.110407239819004,.133031674208145,.146606334841629,.151131221719457,.146606334841629,.133031674208145,.110407239819004,.078733031674208,.038009049773756,-.011764705882353,-.070588235294118]];function xn(n,t,e,i=1){const r=Math.max(1e-6,e*.5),s=Math.max(1e-6,t*i),o=Math.exp(-s*r),a=1-o;if(!Number.isFinite(a)||Math.abs(a)<=1e-9)return Math.abs(n)<=r?1:0;if(Math.abs(n)<r){const l=2-2*o*Math.cosh(s*n),u=2*Math.sinh(s*r)*a;return!Number.isFinite(u)||Math.abs(u)<=1e-9?0:l/u}return Math.exp(-s*Math.abs(n))/a}function Ms(n,t,e,i,r,s){const o=n.length;if(o===0)return[];if(s<1)return n;const a=Math.min(s,32),l=new Array(o).fill(0);l[0]=n[0];for(let d=1;d<o;d++)l[d]=l[d-1]+n[d];const u=(d,f)=>{const p=Math.max(0,d),m=Math.min(o-1,f);return m<p?n[Math.max(0,Math.min(o-1,d))]??0:(l[m]-(p>0?l[p-1]:0))/(m-p+1)},h=a*2,c=a,g=1;for(let d=Math.max(t+c,i-h);d<i;d++){const f=Math.max(g,Math.trunc((i-d)*c/Math.max(1,h)));n[d]=u(d-f,d+f)}for(let d=Math.min(r+h-1,e-c-1);d>r;d--){const f=Math.max(g,Math.trunc((d-r)*c/Math.max(1,h)));n[d]=u(d-f,d+f)}for(let d=c+1;d<i-h;d++)n[d]=u(d-c,d+c);for(let d=Math.min(r+h,e-c-1);d<o-c-1;d++)n[d]=u(d-c,d+c);return n}function ws(n){return!Number.isFinite(n)||Math.abs(n)<=1e-9?1:Math.sin(n)/n}let Hn=null,Di=null;function mo(){if(Hn)return Hn;const n=.625,t=1/128,e=Math.max(16,Math.round(n*2/t)+1),i=[],r=[];for(let c=0;c<e;c++){const g=-n+c*t;i.push(g),r.push(Math.abs(g)<=n?xn(g,yn,.125,1):0)}const s=4,o=1/1024,a=Math.round(s/o)+1,l=new Array(a).fill(0),u=new Array(a).fill(1);let h=0;for(let c=0;c<r.length;c++)h+=r[c];h=Math.max(1e-9,h);for(let c=0;c<a;c++){const g=c*o;l[c]=g;let d=0;for(let f=0;f<i.length;f++)d+=r[f]*Math.cos(2*Math.PI*g*i[f]);u[c]=Math.max(1e-6,Math.abs(d)/h)}return Hn={freqs:l,values:u},Hn}function go(n,t){const e=Math.max(1e-6,ws(Math.PI*n*t)),i=mo(),r=nn(Math.max(0,Math.min(i.freqs[i.freqs.length-1],n)),i.freqs,i.values);return Math.max(1e-6,e*r)}function yo(){if(Di)return Di;const n=new Array(ae/16*4).fill(1),t=ae*16,e=new Float32Array(t);for(let s=0;s<t;s++){const o=(s-t/2)/(16*Et);e[s]=Math.abs(o)<=.625?xn(o,yn,fi,1):0}const i=new In(t);i.transform(e);const r=Math.max(1e-9,Math.abs(i._real[0]));n[0]=1;for(let s=1;s<n.length;s++){const o=ws(Math.PI*s/256),a=Math.max(1e-6,Math.hypot(i._real[s],i._imag[s])/r);n[s]=Math.max(1e-6,o*a)}return Di=n,n}function xo(n,t,e){if(n.length===0||t.length!==n.length||!(e>0))return null;const i=n[0],r=n[n.length-1],s=Math.floor((r-i)/e);if(s<16)return null;const o=Math.max(0,Math.min(s-1,Math.round(-i/e))),a=Math.max(2,Math.round(2/e)),l=Math.max(1,Math.round(.5/e)),u=5,h=.125/Math.max(e,1e-6),c=D=>i+D*e,d=((D,Z)=>{let H=Math.max(0,D),T=Math.min(s-1,Z);if(T-H<8)return null;let Q=0;for(;;){const G=new Array(s).fill(0),Y=new Array(s).fill(0);for(let $=0;$<n.length;$++){const Ft=Math.max(0,Math.min(s-1,Math.trunc((n[$]-i)/e))),Tt=Math.max(H,Ft-u),xt=Math.min(T-1,Ft+u);for(let wt=Tt;wt<=xt;wt++){const ne=c(wt),Jt=Math.max(0,1-Math.abs((n[$]-ne)*1.75*h));Jt<=0||(G[wt]+=t[$]*Jt,Y[wt]+=Jt)}}const et=new Array(s).fill(te);let nt=0,it=0,ct=0,dt=0,rt=-1,Bt=-1;const Xt=Math.max(o-Math.round(s/8),H+a),j=Math.min(o+Math.round(s/8),T-a);for(let $=Math.max(0,H-1);$<=Math.min(s-1,T+1);$++)Y[$]>0&&(et[$]=G[$]/Y[$],$<Xt&&(nt+=et[$],ct++),$>j&&(it+=et[$],dt++),rt<0&&(rt=$),Bt=$);if(rt<0||Bt<0||ct===0||dt===0)return null;for(let $=rt-1;$>=0;$--)et[$]=et[rt];for(let $=Bt+1;$<s;$++)et[$]=et[Bt];const ut=Math.max(2,a),ot=new Array(s).fill(0);let pt=o;for(let $=ut+1;$<s-1-ut;$++){let Ft=0,Tt=0;for(let xt=-ut;xt<=ut;xt++)Tt+=et[$+xt]*xt,Ft+=xt*xt;ot[$]=Ft>0?Tt/Ft:0,Math.abs(ot[$])>Math.abs(ot[pt]??0)&&$>H+ut&&$<T-ut-1&&(pt=$)}const mt=Math.max(1,Math.round(2/e)),ee=Math.max(mt+1,Math.round(12/e)),Dt=Math.abs(pt-o);if(Dt>mt&&Dt<ee)return null;let _t=0;for(let $=Math.max(0,o-ut);$<=Math.min(s-1,o+ut);$++)Math.abs(ot[$])>Math.abs(_t)&&(_t=ot[$]);if(!Number.isFinite(_t)||Math.abs(_t)<=1e-9)return null;const ce=Math.abs(_t*.001);let $t=!1,Ae=!1;for(let $=Math.max(0,o-a);$>=H+l;$--)if(ot[$]*_t<0&&Math.abs(ot[$])>ce){let Ft=0,Tt=0,xt=0;for(let wt=$;wt>=H;wt--)ot[wt]*_t<0&&(Ft++,Tt=Math.max(Tt,Math.abs(ot[wt]))),xt++;if(Ft>xt*.4&&Tt/Math.abs(_t)>.25||Ft>.9*xt&&xt>a){H=Math.min($,o-a),$t=!0;break}}for(let $=Math.min(s-1,o+a);$<T-l;$++)if(ot[$]*_t<0&&Math.abs(ot[$])>ce){let Ft=0,Tt=0,xt=0;for(let wt=$;wt<T;wt++)ot[wt]*_t<0&&(Ft++,Tt=Math.max(Tt,Math.abs(ot[wt]))),xt++;if(Ft>xt*.4&&Tt/Math.abs(_t)>.25||Ft>.9*xt&&xt>a){T=Math.max($,o+a),$t=!0;break}}if($t&&((o-H<Math.max(1,Math.round(4/e))||T-o<Math.max(1,Math.round(4/e)))&&(Ae=!0),Q<2)){Q++;continue}return{sampled:et,fftLeft:H,fftRight:T,leftMean:nt/ct,rightMean:it/dt,peakSlopeIdx:pt,slopes:ot,clipped:$t,dodgy:Ae}}})(0,s-1);if(!d)return null;const f=2,p=Math.max(d.leftMean,d.rightMean),m=Math.min(d.leftMean,d.rightMean);let y=d.fftLeft,x=d.fftLeft,b=1/0,_=1/0;for(let D=d.fftLeft;D<=d.fftRight;D++){let Z=0,H=0;for(let Y=-f;Y<=f;Y++){const et=d.sampled[Math.max(0,Math.min(s-1,D+Y))];Z+=et,H++}const T=Z/Math.max(1,H),Q=Math.abs(T-m-.1*(p-m)),G=Math.abs(T-m-.9*(p-m));Q<b&&(b=Q,y=D),G<_&&(_=G,x=D)}const k=Math.max(4,Math.abs(y-x)*e);if(y<x){const D=y;y=x,x=D}const v=Math.max(l,l+2*Math.round(k/Math.max(e,1e-6)));y+=v,x-=v;const w=Math.max(Math.abs(y-o),Math.abs(x-o),Math.max(a,Math.round(4/Math.max(e,1e-6)))),P=1.85,A=.5,F=new Array(s).fill(0),C=new Array(s).fill(0);for(let D=0;D<n.length;D++){const Z=Math.max(0,Math.min(s-1,Math.trunc((n[D]-i)/e)));let H=5;Math.abs(Z-o)>A*w&&(H=Math.abs(Z-o)>2*A*w?12:7);const T=Math.max(d.fftLeft,Z-H),Q=Math.min(d.fftRight-1,Z+H);if(Q<o-P*w||T>o+P*w){for(let G=T;G<=Q;G++)F[G]+=t[D],C[G]+=1;continue}for(let G=T;G<=Q;G++){let Y=1;if(Math.abs(G-o)<P*w){const et=c(G);if(Math.abs(G-o)<w*A)Y=xn(n[D]-et,yn,e,1);else{const nt=(Math.abs(G-o)/Math.max(1e-6,w)-A)/Math.max(1e-6,P-A),it=1*(1-nt)+.01*nt;Y=xn(n[D]-et,yn,e,it)}}!(Y>0)||!Number.isFinite(Y)||(F[G]+=t[D]*Y,C[G]+=Y)}}const M=new Array(s).fill(0);let S=-1,I=-1;for(let D=0;D<s;D++)C[D]>0?(M[D]=F[D]/C[D],S<0&&(S=D),I=D):M[D]=te;if(S<0||I<0)return null;const N=Math.max(1,Math.round(3/Math.max(e,1e-6)));let E=M[S],X=1;for(let D=S+1;D<o&&X<N;D++)M[D]!==te&&(E+=M[D],X++);E/=X;let U=M[I],L=1;for(let D=I-1;D>o&&L<N;D--)M[D]!==te&&(U+=M[D],L++);U/=L;for(let D=S-1;D>=0;D--)M[D]=E;for(let D=I+1;D<s;D++)M[D]=U;for(let D=S+1;D<I;D++){if(M[D]!==te)continue;let Z=D-1;for(;Z>=0&&M[Z]===te;)Z--;let H=D+1;for(;H<s&&M[H]===te;)H++;if(Z>=0&&H<s){const T=(D-Z)/Math.max(1,H-Z);M[D]=M[Z]*(1-T)+M[H]*T}else Z>=0?M[D]=M[Z]:H<s&&(M[D]=M[H])}const O=[...E<=U?M:[...M].reverse()],z=E<=U?O:O.reverse(),J=Math.max(Math.round(o-P*w),d.fftLeft+2),V=Math.min(Math.round(o+P*w),d.fftRight-3),q=Math.max(1,Math.round(2/Math.max(e,1e-6))),W=Ms(z,d.fftLeft,d.fftRight,J,V,q),tt=new Array(s).fill(0);let at=W[Math.max(0,Math.min(s-1,d.fftLeft))]??W[0]??0;for(let D=d.fftLeft;D<=d.fftRight;D++){const Z=W[D]??at,H=W[Math.min(s-1,D+1)]??Z;tt[D]=H-at,at=Z}return{esf:W,lsfFull:tt}}function Ss(n,t,e=Lt){if(n.length===0||t.length!==n.length)return null;const i=ae,r=i/2,s=fi,o=2*Et,a=Math.max(1,Math.round(.5*Et)),l=5,u=Math.max(0,Math.round(r-e*Et)),h=Math.min(i-1,Math.round(r+e*Et));if(h-u<32)return null;let c=new Array(i).fill(te),g=0,d=0,f=0,p=0,m=u,y=h,x=u,b=h,_=0;for(;;){const T=new Array(i).fill(0),Q=new Array(i).fill(0);c=new Array(i).fill(te),g=0,d=0,f=0,p=0;let G=-1,Y=-1;for(let j=0;j<n.length;j++){const ut=Math.trunc(n[j]*Et+r),ot=Math.max(x,ut-l),pt=Math.min(b-1,ut+l);for(let mt=ot;mt<=pt;mt++){const ee=(mt-r)*s,Dt=Math.max(0,1-Math.abs((n[j]-ee)*1.75));Dt>0&&(Q[mt]+=t[j]*Dt,T[mt]+=Dt)}}const et=Math.max(r-i/8,x+2*Et),nt=Math.min(r+i/8,b-2*Et);for(let j=Math.max(0,x-1);j<=Math.min(i-1,b+1);j++)T[j]>0&&(c[j]=Q[j]/T[j],j<et&&(g+=c[j],f++),j>nt&&(d+=c[j],p++),G<0&&(G=j),Y=j);if(G<0||Y<0||f<=0||p<=0)return null;for(let j=G-1;j>=0;j--)c[j]=c[G];for(let j=Y+1;j<i;j++)c[j]=c[Y];const it=new Array(i).fill(0);let ct=r;const dt=2*Et;for(let j=dt+1;j<i-1-dt;j++){let ut=0,ot=0;for(let pt=-dt;pt<=dt;pt++)ot+=c[j+pt]*pt,ut+=pt*pt;it[j]=ut>0?ot/ut:0,Math.abs(it[j])>Math.abs(it[ct]??0)&&j>x+dt&&j<b-dt-1&&(ct=j)}if(Math.abs(ct-r)>2*Et&&Math.abs(ct-r)<12*Et)return null;let rt=0;for(let j=Math.max(0,r-dt);j<=Math.min(i-1,r+dt);j++)Math.abs(it[j])>Math.abs(rt)&&(rt=it[j]);if(!Number.isFinite(rt)||Math.abs(rt)<=1e-9)return null;const Bt=Math.abs(rt*.001);m=x,y=b;let Xt=!1;for(let j=r-o;j>=x+a;j--)if(it[j]*rt<0&&Math.abs(it[j])>Bt){let ut=0,ot=0,pt=0;for(let mt=j;mt>=x;mt--)it[mt]*rt<0&&(ut++,ot=Math.max(ot,Math.abs(it[mt]))),pt++;if(ut>pt*.4&&ot/Math.abs(rt)>.25||ut>.9*pt&&pt>o){m=Math.min(j,r-o),Xt=!0;break}}for(let j=r+o;j<b-a;j++)if(it[j]*rt<0&&Math.abs(it[j])>Bt){let ut=0,ot=0,pt=0;for(let mt=j;mt<b;mt++)it[mt]*rt<0&&(ut++,ot=Math.max(ot,Math.abs(it[mt]))),pt++;if(ut>pt*.4&&ot/Math.abs(rt)>.25||ut>.9*pt&&pt>o){y=Math.max(j,r+o),Xt=!0;break}}if(Xt&&_<2){x=m,b=y,_++;continue}break}const k=Math.max(g/f,d/p),v=Math.min(g/f,d/p);let w=m,P=m,A=1/0,F=1/0;for(let T=m;T<=y;T++){const Q=(c[Math.max(0,T-2)]+c[Math.max(0,T-1)]+c[T]+c[Math.min(i-1,T+1)]+c[Math.min(i-1,T+2)])/5,G=Math.abs(Q-v-.1*(k-v)),Y=Math.abs(Q-v-.9*(k-v));G<A&&(A=G,w=T),Y<F&&(F=Y,P=T)}if(w<P){const T=w;w=P,P=T}const C=Math.max(4,Math.abs(w-P)*s),M=Math.max(a,a+2*Math.trunc(C/Math.max(s,1e-6)));w+=M,P-=M;const S=Math.max(Math.abs(w-r),Math.abs(P-r),Math.max(o,Math.trunc(4/Math.max(s,1e-6)))),I=new Array(i).fill(0),N=new Array(i).fill(0),E=1.85,X=.5;for(let T=0;T<n.length;T++){const Q=Math.trunc(n[T]*Et+r);let G=5;Math.abs(Q-r)>X*S&&(G=Math.abs(Q-r)>2*X*S?12:7);const Y=Math.max(m,Q-G),et=Math.min(y-1,Q+G);if(et<r-E*S||Y>r+E*S){for(let nt=Y;nt<=et;nt++)I[nt]+=t[T],N[nt]+=1;continue}for(let nt=Y;nt<=et;nt++){let it=1;if(Math.abs(nt-r)<E*S){const ct=(nt-r)*s;if(Math.abs(nt-r)<S*X)it=xn(n[T]-ct,yn,s,1);else{const dt=(Math.abs(nt-r)/Math.max(1e-6,S)-X)/Math.max(1e-6,E-X),rt=1*(1-dt)+.01*dt;it=xn(n[T]-ct,yn,s,rt)}}!(it>0)||!Number.isFinite(it)||(I[nt]+=t[T]*it,N[nt]+=it)}}const U=new Array(i).fill(0);let L=-1,B=-1;for(let T=Math.max(0,m-1);T<=Math.min(i-1,y+1);T++)N[T]>0?(U[T]=I[T]/N[T],L<0&&(L=T),B=T):U[T]=te;if(L<0||B<0)return null;const O=3*Et;let z=U[L],J=1;for(let T=L+1;T<r&&J<O;T++)U[T]!==te&&(z+=U[T],J++);z/=Math.max(1,J);let V=U[B],q=1;for(let T=B-1;T>r&&q<O;T--)U[T]!==te&&(V+=U[T],q++);V/=Math.max(1,q);for(let T=L-1;T>=0;T--)U[T]=z;for(let T=B+1;T<i;T++)U[T]=V;const W=Math.max(Math.trunc(r-E*S),m+2),tt=Math.min(Math.trunc(r+E*S),y-3),at=Math.max(1,Math.trunc(2/Math.max(s,1e-6))),D=Ms(U,m,y,W,tt,at),Z=new Array(i).fill(0);let H=D[Math.max(0,Math.min(i-1,m))]??D[0]??0;for(let T=m;T<=y;T++){const Q=D[T]??H;Z[T]=(D[Math.min(i-1,T+1)]??Q)-H,H=Q}return{esf:D,lsfFull:Z}}function bo(n){const t=new Array(n.length).fill(0);if(n.length===0)return t;t[0]=n[0];for(let e=1;e<n.length;e++){const i=rr(n[e]-n[e-1]);t[e]=t[e-1]+i}return t}function rr(n){if(!Number.isFinite(n))return n;let t=(n+Math.PI)%(2*Math.PI);return t<0&&(t+=2*Math.PI),t-Math.PI}function _o(n,t,e=0){if(n.length===0)return[];const i=Number.isFinite(e)?e:0,r=n.map((o,a)=>{const l=t[a]??0,u=-2*Math.PI*i*l;return rr(o-u)});return bo(r).map((o,a)=>{const l=t[a]??0,u=-2*Math.PI*i*l;return o+u})}function Mo(n,t,e,i=.05,r=Number.POSITIVE_INFINITY){const s=Math.min(n.length,t.length);if(s<2)return null;let o=0;if(e)for(let p=1;p<s;p++){const m=e[p];Number.isFinite(m)&&m>o&&(o=m)}const a=e&&o>0?Math.max(1e-6,o*i):0;let l=0,u=0,h=0,c=0,g=0,d=0;for(let p=1;p<s;p++){const m=n[p],y=t[p];if(!Number.isFinite(m)||!Number.isFinite(y)||Math.abs(m)<=1e-12||m>r)continue;const x=e?e[p]:1;if(!Number.isFinite(x)||x<=a)continue;const b=e?x*x:1;l+=b,u+=b*m,h+=b*y,c+=b*m*m,g+=b*m*y,d++}if(d<4||l<=0)return null;const f=l*c-u*u;return Math.abs(f)<=1e-12?null:{slope:(l*g-u*h)/f,intercept:(h*c-u*g)/f,used:d,threshold:a}}function wo(n,t,e=Number.POSITIVE_INFINITY){const i=[],r=[],s=Math.min(n.length,t.length);for(let o=1;o<s;o++)Number.isFinite(n[o])&&Number.isFinite(t[o])&&Math.abs(n[o])>1e-12&&n[o]<=e&&(i.push(n[o]),r.push(t[o]));return i.length<2?{slope:0,intercept:Number.isFinite(t[0])?t[0]:0,used:i.length}:{...ir(i,r),used:i.length}}function So(n,t,e,i=.05,r=Number.POSITIVE_INFINITY,s=0){const o=Math.min(n.length,t.length);if(o<4)return null;let a=0;if(e)for(let b=1;b<o;b++){const _=e[b];Number.isFinite(_)&&_>a&&(a=_)}const l=e&&a>0?Math.max(1e-6,a*i):0,u=[];for(let b=1;b<o;b++){const _=n[b],k=t[b];if(!Number.isFinite(_)||!Number.isFinite(k)||Math.abs(_)<=1e-12||_>r)continue;const v=e?e[b]:1;!Number.isFinite(v)||v<=l||u.push({freq:_,phase:k,weight:e?v*v:1})}if(u.length<4)return null;const h=b=>{let _=0,k=0;for(const F of u){const C=F.phase+2*Math.PI*b*F.freq;_+=F.weight*Math.sin(C),k+=F.weight*Math.cos(C)}const v=Math.atan2(_,k);let w=0,P=0;const A=.65;for(const F of u){const C=F.phase+2*Math.PI*b*F.freq,M=Math.abs(rr(C-v)),S=M<=A?M*M:A*(2*M-A);w+=F.weight*S,P+=F.weight}return{score:P>0?w/P:Number.POSITIVE_INFINITY,intercept:v}},c=Number.isFinite(s)?s:0,g=Math.max(2,Math.min(8,Math.abs(c)>1e-6?4:2)),d=.02;let f=c,p=h(f);for(let b=c-g;b<=c+g+d*.5;b+=d){const _=h(b);_.score<p.score&&(p=_,f=b)}let m=f-d*2,y=f+d*2;for(let b=0;b<32;b++){const _=m+(y-m)/3,k=y-(y-m)/3,v=h(_).score,w=h(k).score;v<w?y=k:m=_}const x=(m+y)*.5;return p=h(x),{slope:-2*Math.PI*x,intercept:p.intercept,used:u.length,threshold:l}}function vo(n,t){if(Number.isFinite(t)&&t>0)return t;let e=0;for(const i of n)Number.isFinite(i)&&i>e&&(e=i);return e>0?e:Number.POSITIVE_INFINITY}function vs(n,t,e,i=Number.POSITIVE_INFINITY,r=0){const s=_o(n,t,r),o=vo(t,i),a=So(t,n,e,.05,o,r),l=a?null:Mo(t,s,e,.05,o),u=a||l?null:wo(t,s,o),h=(a==null?void 0:a.slope)??(l==null?void 0:l.slope)??(u==null?void 0:u.slope)??0,c=(a==null?void 0:a.intercept)??(l==null?void 0:l.intercept)??(u==null?void 0:u.intercept)??0,g=s.map((x,b)=>x-(h*(t[b]??0)+c)),d=Number.isFinite(g[0])?g[0]:0,f=t.map(x=>h*x+c+d),p=g.map(x=>x-d),m=c+d,y=Number.isFinite(h)?-h/(2*Math.PI):null;return{raw:[...n],unwrapped:s,linear:f,residual:p,fit:{groupDelayPx:y===null?null:y-r,absoluteGroupDelayPx:y,referenceDelayPx:r,slopeRadPerCycle:Number.isFinite(h)?h:null,interceptRad:Number.isFinite(m)?m:null,fitPointCount:(a==null?void 0:a.used)??(l==null?void 0:l.used)??(u==null?void 0:u.used)??0,fitWeightThreshold:(a==null?void 0:a.threshold)??(l==null?void 0:l.threshold)??0,fitDomain:"cycles-per-pixel",fitMaxFreqCyclesPerPixel:o}}}function Ps(n,t,e){const i=[],r=[],s=[],o=[];for(let a=0;a<e.length;a++){const l=e[a];i.push(nn(l,t,n.raw)),r.push(nn(l,t,n.unwrapped)),s.push(nn(l,t,n.linear)),o.push(nn(l,t,n.residual))}return{ptfRaw:i,ptfUnwrapped:r,ptfLinear:s,ptfResidual:o}}function Cs(n,t){const e=n.map((a,l)=>({dist:a,value:t[l]})).filter(a=>Number.isFinite(a.dist)&&Number.isFinite(a.value)).sort((a,l)=>a.dist-l.dist);if(e.length===0)return{dists:[],vals:[]};const i=Math.max(1,Math.min(16,Math.floor(e.length*.1)));let r=0,s=0;for(let a=0;a<i;a++)r+=e[a].value,s+=e[e.length-1-a].value;r/=i,s/=i;const o=r<=s?e:e.map(a=>({dist:-a.dist,value:a.value})).sort((a,l)=>a.dist-l.dist);return{dists:o.map(a=>a.dist),vals:o.map(a=>a.value)}}function Po(n,t,e){const i=dn(e,"RAW green-site sampling"),r=t%2,s=n%2;return i==="RGGB"||i==="BGGR"?(r+s)%2!==0:i==="GBRG"||i==="GRBG"?(r+s)%2===0:!1}function Co(n,t){return n+t&1?2:1}function Gr(n,t){return(t&1)<<1|n&1}function Mn(n,t,e){return n===void 0?0:typeof n=="number"?Number.isFinite(n)?n:0:Number.isFinite(n[Gr(t,e)])?n[Gr(t,e)]:0}function vt(n,t,e,i){return i!==void 0&&i!=="default"?Co(n,t)===i:Po(n,t,e)}function ri(n){return n.length===0?0:n.reduce((t,e)=>t+e,0)/n.length}function ko(n,t,e){const i=(e%t+t)%t,r=Math.floor(i),s=(r+1)%t,o=i-r,a=r<n.length?n[r]:0,l=s<n.length?n[s]:0;return a*(1-o)+l*o}function ks(n,t){const e=n.length;if(e===0)return 0;if(t<=0)return n[0];if(t>=e-1)return n[e-1];const i=Math.floor(t),r=Math.min(e-1,i+1),s=t-i;return n[i]*(1-s)+n[r]*s}function Fo(n,t,e){const i=n.length,r=new Array(i).fill(0);for(let s=0;s<i;s++)r[s]=ks(n,s-e+t);return r}function si(n,t,e){const i=n.length;if(i===0)return{peakPos:0,peakIdx:0,peakVal:0};const r=Math.max(0,Math.floor(t-e)),s=Math.min(i-1,Math.ceil(t+e));let o=Math.max(0,Math.min(i-1,Math.round(t))),a=-1/0;for(let u=r;u<=s;u++){const h=Math.abs(n[u]);h>a&&(a=h,o=u)}Number.isFinite(a)||(a=Math.abs(n[o]??0));let l=o;if(o>0&&o<i-1){const u=n[o]>=0?1:-1,h=u*n[o-1],c=u*n[o],g=u*n[o+1],d=h-2*c+g;if(Number.isFinite(d)&&Math.abs(d)>1e-9){const f=.5*(h-g)/d;Number.isFinite(f)&&Math.abs(f)<=1&&(l=o+f)}}return{peakPos:l,peakIdx:o,peakVal:Math.abs(ks(n,l))}}function Ao(n,t,e,i,r,s,o,a,l){const u=Math.floor(i.x),h=Math.floor(i.y),c=Math.floor(i.w),g=Math.floor(i.h),d=[],f=(p,m)=>{if(p<0||m<0||p>=t||m>=e)return null;const y=r+p,x=s+m;return Math.max(0,n[m*t+p]-Mn(l,y,x))};for(let p=0;p<g;p++){const m=[],y=h+p;for(let x=0;x<c;x++){const b=u+x,_=r+b,k=s+y;if(vt(_,k,o,a)){m.push(f(b,y)??0);continue}const v=[],w=f(b-1,y),P=f(b+1,y),A=f(b,y-1),F=f(b,y+1);if(w!==null&&vt(_-1,k,o,a)&&v.push(w),P!==null&&vt(_+1,k,o,a)&&v.push(P),A!==null&&vt(_,k-1,o,a)&&v.push(A),F!==null&&vt(_,k+1,o,a)&&v.push(F),v.length===0){const C=[],M=f(b-1,y-1),S=f(b+1,y-1),I=f(b-1,y+1),N=f(b+1,y+1);M!==null&&vt(_-1,k-1,o,a)&&C.push(M),S!==null&&vt(_+1,k-1,o,a)&&C.push(S),I!==null&&vt(_-1,k+1,o,a)&&C.push(I),N!==null&&vt(_+1,k+1,o,a)&&C.push(N),m.push(ri(C));continue}m.push(ri(v))}d.push(m)}return d}function Fs(n,t,e,i,r,s,o,a,l){const u=Math.floor(i.x),h=Math.floor(i.y),c=Math.floor(i.w),g=Math.floor(i.h),d=[];for(let f=0;f<g;f++){const p=h+f;for(let m=0;m<c;m++){const y=u+m,x=r+y,b=s+p;vt(x,b,o,a)&&d.push({x,y:b,value:Math.max(0,n[p*t+y]-Mn(l,x,b))})}}return d}function As(n,t,e,i,r){const s=Math.floor(i.x),o=Math.floor(i.y),a=Math.floor(i.w),l=Math.floor(i.h),u=(r==null?void 0:r.globalX)??0,h=(r==null?void 0:r.globalY)??0,c=!!(r!=null&&r.isThreePlane)&&n.length>=t*e*3,g=r==null?void 0:r.threePlaneChannel,d=[];for(let f=0;f<l;f++){const p=o+f,m=p*t;for(let y=0;y<a;y++){const x=s+y;let b=0;if(!c)b=Math.max(0,n[m+x]-Mn(r==null?void 0:r.blackLevel,u+x,h+p));else{const _=(m+x)*3;if(g!==void 0)b=n[_+g];else{const k=n[_],v=n[_+1],w=n[_+2];b=.2126*k+.7152*v+.0722*w}}d.push({x:u+x,y:h+p,value:b})}}return d}function To(n){var s;const t=n.length,e=((s=n[0])==null?void 0:s.length)??0;let i=0,r=0;for(let o=1;o<t-1;o++)for(let a=1;a<e-1;a++)i+=Math.abs(n[o][a+1]-n[o][a-1]),r+=Math.abs(n[o+1][a]-n[o-1][a]);return{gx:i,gy:r}}function Ts(n,t,e,i,r,s,o){var g;const a=n.length,l=((g=n[0])==null?void 0:g.length)??0,u=(d,f,p,m,y,x,b)=>{const _=b?a:l,k=Math.max(0,f-3),v=Math.min(_,f+4);let w=0,P=0;for(let F=k;F<v;F++)w+=d[F],P+=F*d[F];if(w<=0)return null;const A=P/w;return b?{x:t+m*y,y:p+A*x,weight:w}:{x:p+A*x,y:e+m*y,weight:w}},h=(d,f,p,m,y,x,b)=>{const _=Math.max(3,Math.min(Math.max(3,Math.floor(p/3)),Math.max(4,Math.round(p*.12)))),k=d.map(S=>{let I=-1/0,N=-1;for(let E=0;E<S.length;E++)S[E]>I&&(I=S[E],N=E);return{peakValue:I,peakIndex:N}}),v=(f-1)*.5,w=k.map((S,I)=>({...S,index:I})).filter(S=>S.peakValue>1&&S.peakIndex>=0).sort((S,I)=>{const N=I.peakValue-S.peakValue;return Math.abs(N)>1e-6?N:Math.abs(S.index-v)-Math.abs(I.index-v)});if(w.length===0)return[];const P=w[0],A=new Array(f).fill(null),F=u(d[P.index],P.peakIndex,m,P.index,y,x,b);if(!F)return[];A[P.index]=F;const C=(S,I)=>{const N=d[S],E=k[S];if(!(E.peakValue>1)||E.peakIndex<0)return null;const X=Math.max(0,Math.floor(I-_)),U=Math.min(N.length,Math.ceil(I+_+1));let L=-1/0,B=-1;for(let V=X;V<U;V++)N[V]>L&&(L=N[V],B=V);if(B<0||!(L>1))return null;const O=Math.max(1e-6,E.peakValue),z=Math.abs(B-I)<=_,J=L>=O*.25;return!z||!J?null:u(N,B,m,S,y,x,b)};let M=F?b?(F.y-m)/x:(F.x-m)/x:P.peakIndex;for(let S=P.index+1;S<f;S++){const I=C(S,M);I&&(A[S]=I,M=b?(I.y-m)/x:(I.x-m)/x)}M=F?b?(F.y-m)/x:(F.x-m)/x:P.peakIndex;for(let S=P.index-1;S>=0;S--){const I=C(S,M);I&&(A[S]=I,M=b?(I.y-m)/x:(I.x-m)/x)}return A.filter(S=>!!S)};if(s){const d=n.map(f=>f.map((p,m)=>m===0?0:Math.abs(p-f[m-1])));return h(d,a,l,t,r,i,!1)}const c=[];for(let d=0;d<l;d++){const f=new Array(a).fill(0);for(let p=1;p<a;p++)f[p]=Math.abs(n[p][d]-n[p-1][d]);c.push(f)}return h(c,l,a,e,i,r,!0)}function Ce(n){if(n.length<2)return null;let t=0,e=0,i=0;for(const c of n)t+=c.weight,e+=c.x*c.weight,i+=c.y*c.weight;if(t<=0)return null;e/=t,i/=t;let r=0,s=0,o=0;for(const c of n){const g=c.x-e,d=c.y-i;r+=c.weight*g*g,s+=c.weight*d*d,o+=c.weight*g*d}r/=t,s/=t,o/=t;const a=.5*Math.atan2(2*o,r-s);let l=Math.cos(a),u=Math.sin(a);const h=Math.hypot(l,u);return!Number.isFinite(h)||h<=1e-9?null:(l/=h,u/=h,(l<0||Math.abs(l)<=1e-9&&u<0)&&(l=-l,u=-u),{pointX:e,pointY:i,dirX:l,dirY:u,orientation:Math.abs(l)>=Math.abs(u)?1:2})}function Io(n,t){if(n.length!==4||t.length!==4||n.some(i=>i.length!==4))return null;const e=n.map((i,r)=>[...i,t[r]]);for(let i=0;i<4;i++){let r=i,s=Math.abs(e[i][i]);for(let a=i+1;a<4;a++){const l=Math.abs(e[a][i]);l>s&&(s=l,r=a)}if(!(s>1e-12))return null;if(r!==i){const a=e[i];e[i]=e[r],e[r]=a}const o=e[i][i];for(let a=i;a<=4;a++)e[i][a]/=o;for(let a=0;a<4;a++){if(a===i)continue;const l=e[a][i];if(!(Math.abs(l)<=1e-12))for(let u=i;u<=4;u++)e[a][u]-=l*e[i][u]}}return[e[0][4],e[1][4],e[2][4],e[3][4]]}function Ro(n){if(n.length<4)return 0;const t=[...n].sort((g,d)=>g.x-d.x),e=t[0].x,r=t[t.length-1].x-e;if(!(r>1e-6))return 0;const s=16,o=[];for(let g=0;g<s;g++){const d=Math.max(0,Math.floor((g-1.5)*t.length/s)),f=Math.min(t.length-1,Math.floor((g+2.5)*t.length/s));if(f<d)continue;let p=0,m=0,y=0;for(let x=d;x<=f;x++)p+=t[x].x,m+=t[x].y,y++;y>0&&o.push({x:p/y,y:m/y})}if(o.length<4)return 0;const a=[.05952381,0,-.03571429,-.04761905,-.03571429,0,.05952381],l=new Array(o.length).fill(0),u=3;for(let g=0;g<o.length;g++){let d=0;for(let f=-u;f<=u;f++){const p=g+f;p<0||p>=o.length||(d+=a[f+u]*o[p].y)}l[g]=d}let h=0,c=1/0;for(let g=0;g<l.length-1;g++){const d=l[g],f=l[g+1];if(d===0){const b=Math.abs(o[g].x);b<c&&(c=b,h=o[g].x);continue}if(d*f>=0)continue;const p=f-d,m=Math.abs(p)>1e-12?-d/p:.5,y=o[g].x+(o[g+1].x-o[g].x)*m,x=Math.abs(y);x<c&&(c=x,h=y)}return!Number.isFinite(h)||h<e+.3*r||h>e+.7*r?0:h}function No(n){if(n.length<8)return null;const t=[...n].filter(f=>Number.isFinite(f.x)&&Number.isFinite(f.y)&&Number.isFinite(f.weight)).sort((f,p)=>f.x-p.x);if(t.length<8)return null;const i=[Ro(t),0,.5*(t[Math.floor((t.length-1)*.5)].x+t[Math.ceil((t.length-1)*.5)].x)];let r=null;for(const f of i){if(!Number.isFinite(f))continue;let p=0,m=0;for(const y of t)y.x<=f?p++:m++;if(p>=4&&m>=4){r=f;break}}if(r===null)return null;const s=Array.from({length:4},()=>new Array(4).fill(0)),o=new Array(4).fill(0);for(const f of t){const p=f.x,m=f.y,y=Math.max(1e-6,f.weight),x=p<=r?[p*p,p,1,0]:[2*r*p-r*r,p,1,(p-r)*(p-r)];for(let b=0;b<4;b++){o[b]+=y*x[b]*m;for(let _=0;_<4;_++)s[b][_]+=y*x[b]*x[_]}}const a=Io(s,o);if(!a)return null;const[l,u,h,c]=a,g=u+2*(l-c)*r,d=h+(c-l)*r*r;return[l,u,h,c,g,d].every(f=>Number.isFinite(f))?{splitX:r,left:[l,u,h],right:[c,g,d]}:null}function Lo(n,t,e){const[i,r,s]=e;if(Math.abs(i)<=1e-12){const b=1+r*r;return b>1e-12?[(n-r*(s-t))/b]:[n]}const o=2*i*i,a=3*i*r,l=1+2*i*s-2*i*t+r*r,u=r*s-t*r-n;if(Math.abs(o)<=1e-12)return[n];const h=a/o,c=l/o,g=u/o,d=(h*h-3*c)/9,f=(2*h*h*h-9*h*c+27*g)/54,p=f*f-d*d*d;if(p<0&&d>0){const b=Math.acos(Math.max(-1,Math.min(1,f/Math.sqrt(d*d*d)))),_=-2*Math.sqrt(d);return[_*Math.cos(b/3)-h/3,_*Math.cos((b+2*Math.PI)/3)-h/3,_*Math.cos((b-2*Math.PI)/3)-h/3]}const m=Math.sqrt(Math.max(0,p)),y=-Math.sign(f||1)*Math.cbrt(Math.abs(f)+m),x=Math.abs(y)<=1e-12?0:d/y;return[y+x-h/3]}function Is(n,t){if(n.length<8)return null;const e=-t.dirY,i=t.dirX,r=n.map(o=>({x:(o.x-t.pointX)*t.dirX+(o.y-t.pointY)*t.dirY,y:(o.x-t.pointX)*e+(o.y-t.pointY)*i,weight:o.weight})),s=No(r);return s?{...t,normalX:e,normalY:i,splitX:s.splitX,left:s.left,right:s.right}:null}function Eo(n,t){const e=n.x-t.pointX,i=n.y-t.pointY,r=e*t.dirX+i*t.dirY,s=e*t.normalX+i*t.normalY,o=r<t.splitX?t.left:t.right,a=Lo(r,s,o);let l=s,u=Number.POSITIVE_INFINITY;for(const h of a){if(!Number.isFinite(h))continue;const c=o[0]*h*h+o[1]*h+o[2],g=r-h,d=s-c,f=Math.hypot(g,d);Number.isFinite(f)&&f<u&&(u=f,l=(d>=0?1:-1)*f)}return Number.isFinite(u)?l:s}function Rs(n,t,e,i,r,s,o,a=Lt){if(!t||t.length<8||n.length===0)return null;const l=t.filter(S=>Number.isFinite(S.x)&&Number.isFinite(S.y)).map(S=>({x:S.x,y:S.y,weight:1}));if(l.length<8)return null;const u=Ce(l);if(!u)return null;const h=e.p2.x-e.p1.x,c=e.p2.y-e.p1.y,g=Math.hypot(h,c);if(!Number.isFinite(g)||g<=1e-6)return null;let d=u.dirX,f=u.dirY;d*h+f*c<0&&(d=-d,f=-f);const p={...u,dirX:d,dirY:f},m=Is(l,p),y=h/g,x=c/g,b=-x,_=y,k=(e.p1.x+e.p2.x)*.5,v=(e.p1.y+e.p2.y)*.5,w=-p.dirY,P=p.dirX,A=Math.abs(y)>=Math.abs(x)?1:2,F=[],C=[];for(const S of n){const I=S.x-k,N=S.y-v,E=I*y+N*x;if(Math.abs(E)>i)continue;const X=I*b+N*_;if(Math.abs(X)>r)continue;const U=m?Eo(S,m):(S.x-p.pointX)*w+(S.y-p.pointY)*P;Number.isFinite(U)&&(F.push(U),C.push(S.value))}if(F.length<8)return null;const M=o?s!=null&&s.forceLegacyModel?Rn(F,C,A,a,r*2):wn(F,C,A,a):pi(F,C,Math.max(2,r*2),s==null?void 0:s.manualBinSize,A,s==null?void 0:s.preferAutoPerEdgeBin);return M?(M.quadraticProjectionUsed=!!m,M):null}function zr(n){if(n.length<2)return null;const t=n.filter(e=>Number.isFinite(e.x)&&Number.isFinite(e.y)).map(e=>({x:e.x,y:e.y,weight:1}));return t.length<2?null:Fn(t,Ce(t))}function Fn(n,t){if(!t||n.length<2)return null;let e=1/0,i=-1/0;for(const s of n){const o=(s.x-t.pointX)*t.dirX+(s.y-t.pointY)*t.dirY;e=Math.min(e,o),i=Math.max(i,o)}if(!Number.isFinite(e)||!Number.isFinite(i))return null;const r=Math.max(.5,(i-e)*.03);return{p1:{x:t.pointX+t.dirX*(e-r),y:t.pointY+t.dirY*(e-r)},p2:{x:t.pointX+t.dirX*(i+r),y:t.pointY+t.dirY*(i+r)}}}function Uo(n,t,e,i,r,s,o,a){return[Fn(n,t),Fn(e,i),Fn(r,s),Fn(o,a)]}function Bo(n,t,e){if(!n||n.length<8)return;const i=n.filter(x=>Number.isFinite(x.x)&&Number.isFinite(x.y)).map(x=>({x:x.x,y:x.y,weight:1}));if(i.length<8)return;const r=Ce(i);if(!r)return;const s=t.p2.x-t.p1.x,o=t.p2.y-t.p1.y,a=Math.hypot(s,o);if(!Number.isFinite(a)||a<=1e-6)return;let l=r.dirX,u=r.dirY;l*s+u*o<0&&(l=-l,u=-u);const h=Is(i,{...r,dirX:l,dirY:u});if(!h)return;const c=n.map(x=>(x.x-h.pointX)*h.dirX+(x.y-h.pointY)*h.dirY).filter(x=>Number.isFinite(x)),g=(t.p1.x-h.pointX)*h.dirX+(t.p1.y-h.pointY)*h.dirY,d=(t.p2.x-h.pointX)*h.dirX+(t.p2.y-h.pointY)*h.dirY;if(Number.isFinite(g)&&c.push(g),Number.isFinite(d)&&c.push(d),c.length<2)return;const f=Math.min(...c),p=Math.max(...c);if(!Number.isFinite(f)||!Number.isFinite(p)||p-f<=1e-6)return;const m=Math.max(21,e),y=[];for(let x=0;x<m;x++){const b=m===1?.5:x/(m-1),_=f+(p-f)*b,k=_<h.splitX?h.left:h.right,v=k[0]*_*_+k[1]*_+k[2];y.push({x:h.pointX+_*h.dirX+v*h.normalX,y:h.pointY+_*h.dirY+v*h.normalY})}return y}function pi(n,t,e,i,r,s=!1,o=!1,a=!1){if(n.length===0||t.length!==n.length)return null;const l=Cs(n,t),u=l.dists,h=l.vals;if(u.length===0)return null;const c=()=>{const w=e/2;let P=0;for(const F of n)Math.abs(F)<=w&&P++;if(P<=0)return .125;const A=40*w/P;return Math.max(.01,Math.min(.125,A))},g=(w,P,A,F,C)=>{if(!(F>0)||!(C>0)||!(A>P))return!1;const M=Math.floor((A-P)/F);if(M<2)return!1;const S=Math.max(P,-C),I=Math.min(A,C);if(!(I>S))return!1;const N=Math.max(0,Math.floor((S-P)/F)),E=Math.min(M,Math.ceil((I-P)/F));if(E<=N)return!1;const X=new Array(E-N).fill(0),U=P+N*F,L=P+E*F;for(let B=0;B<w.length;B++){const O=w[B];if(O<U)continue;if(O>=L)break;const z=Math.floor((O-P)/F);z>=N&&z<E&&X[z-N]++}return X.every(B=>B>0)},d=()=>{const w=u[0],P=u[u.length-1],A=Math.max(0,e*.25),F=.125,C=.5,M=.001,S=Math.round((C-F)/M);for(let I=0;I<=S;I++){const N=Number((F+I*M).toFixed(3));if(g(u,w,P,N,A))return N}return C};let f=.125;i&&i>0?f=Math.max(.01,Math.min(.5,i)):s?f=d():f=c();const p=u[0],m=u[u.length-1],y=Math.floor((m-p)/f);if(y<2)return null;const x=()=>{const w=new Array(y).fill(0),P=new Array(y).fill(0);for(let F=0;F<u.length;F++){const C=(u[F]-p)/f;if(Number.isFinite(C))if(o){const M=Math.floor(C),S=C-M,I=1-S,N=S;M>=0&&M<y&&(w[M]+=h[F]*I,P[M]+=I);const E=M+1;E>=0&&E<y&&(w[E]+=h[F]*N,P[E]+=N)}else{const M=Math.floor(C);M>=0&&M<y&&(w[M]+=h[F],P[M]++)}}let A=h[0];for(let F=0;F<y;F++)P[F]>0?(w[F]/=P[F],A=w[F]):w[F]=A;return w},b=a?null:xo(u,h,f),_=(b==null?void 0:b.esf)??x(),k=(b==null?void 0:b.lsfFull)??fo(_,f),v=Math.max(0,Math.min(y-1,-p/f-.5));return{esf:_,lsfFull:k,binSize:f,orientation:r,zeroIndex:v,shortSidePx:e,fallbackUsed:a||!b}}function wn(n,t,e,i=Lt){if(n.length===0||t.length!==n.length)return null;const r=Cs(n,t),s=r.dists.map((a,l)=>({dist:a,value:r.vals[l]})).filter(a=>Math.abs(a.dist)<i);if(s.length<8)return null;const o=Ss(s.map(a=>a.dist),s.map(a=>a.value),i);return o?{esf:o.esf,lsfFull:o.lsfFull,binSize:fi,orientation:e,zeroIndex:ae/2,shortSidePx:i*2,fallbackUsed:!1,mtfmapperLike:!0,mtfmapperOrderedDists:s.map(a=>a.dist),mtfmapperOrderedVals:s.map(a=>a.value),mtfmapperEffectiveMaxDot:i}:null}function Rn(n,t,e,i=Lt,r=i*2){if(n.length===0||t.length!==n.length)return null;const s=[],o=[];for(let a=0;a<n.length;a++){const l=n[a],u=t[a];!Number.isFinite(l)||!Number.isFinite(u)||Math.abs(l)>=i||(s.push(l),o.push(u))}return s.length<8?null:pi(s,o,Math.max(2,r),void 0,e,!0,!0,!0)}function Pe(n,t,e,i,r,s,o,a,l=Lt){if(s<=0||o<=0)return null;const h=!(!!(a!=null&&a.isThreePlane)&&n.length>=t*e*3)&&((a==null?void 0:a.greenOnly)??!1),c=h?dn(a==null?void 0:a.bayerPattern,"constrained RAW SFR sampling"):null,g=r.p2.x-r.p1.x,d=r.p2.y-r.p1.y,f=Math.hypot(g,d);if(!Number.isFinite(f)||f<=1e-6)return null;const p=g/f,m=d/f,y=-m,x=p,b=(r.p1.x+r.p2.x)*.5,_=(r.p1.y+r.p2.y)*.5,k=Math.abs(p)>=Math.abs(m)?1:2,v=h?Fs(n,t,e,i,0,0,c,a==null?void 0:a.greenPhase,a==null?void 0:a.blackLevel):As(n,t,e,i,{...a,globalX:0,globalY:0});if(v.length===0)return null;if(!(a!=null&&a.disableQuadraticProjection)){const A=Rs(v,a==null?void 0:a.quadraticFitPoints,r,s,o,a,!0,l);if(A)return A}const w=[],P=[];for(const A of v){const F=A.x-b,C=A.y-_,M=F*p+C*m;if(Math.abs(M)>s)continue;const S=F*y+C*x;Math.abs(S)>o||(w.push(S),P.push(A.value))}return w.length<8?null:a!=null&&a.forceLegacyModel?Rn(w,P,k,l):wn(w,P,k,l)}function Do(n,t,e=0){const i=[...n.lsfFull];if(i.length<3)return!1;const r=Math.max(n.binSize,1e-6),s=Number.isFinite(n.zeroIndex)?n.zeroIndex:i.length/2,o=Math.max(1,Math.round((n.shortSidePx??0)*.5/r));let{peakPos:a,peakIdx:l,peakVal:u}=si(i,s,o);const h=u*.2;let c=0,g=i.length-1;for(let p=l;p>=0;p--)if(i[p]<h){c=p;break}for(let p=l;p<i.length;p++)if(i[p]<h){g=p;break}const d=g-c;if(t&&d>0){const p=d*4,m=[],y=[];if(e>0){const x=Math.max(0,l-p-e),b=Math.max(0,l-p);for(let v=x;v<b;v++)m.push(v),y.push(i[v]);const _=Math.min(i.length,l+p),k=Math.min(i.length,l+p+e);for(let v=_;v<k;v++)m.push(v),y.push(i[v])}else{for(let x=0;x<Math.max(0,l-p);x++)m.push(x),y.push(i[x]);for(let x=Math.min(i.length,l+p);x<i.length;x++)m.push(x),y.push(i[x])}if(m.length>2){const{slope:x,intercept:b}=ir(m,y);for(let _=0;_<i.length;_++)i[_]=i[_]-(x*_+b);({peakPos:a}=si(i,s,o))}}return Math.abs(a-s)*r<=Math.max(1e-6,(n.shortSidePx??0)/6)}function Oo(n){const t=n.length;if(t<3)return!1;let e=0,i=-1/0;for(let o=0;o<t;o++){const a=Math.abs(n[o]);a>i&&(i=a,e=o)}const r=t/3,s=2*t/3;return e>=r&&e<=s}function Mt(n,t,e){const i=Math.max(0,Math.floor(n.x)),r=Math.max(0,Math.floor(n.y)),s=Math.min(t,Math.ceil(n.x+n.w)),o=Math.min(e,Math.ceil(n.y+n.h)),a=s-i,l=o-r;return a<2||l<2?null:{x:i,y:r,w:a,h:l}}function sr(n,t,e,i){const r=[],s=n.x,o=n.y,a=n.x+n.w,l=n.y+n.h,u=n.x+n.w*.5,h=n.y+n.h*.5,c=[{x:s,y:o},{x:a,y:o},{x:a,y:l},{x:s,y:l},{x:u,y:o},{x:a,y:h},{x:u,y:l},{x:s,y:h},{x:u,y:h}];for(const g of c){const d=Fe(g,t);Number.isFinite(d.x)&&Number.isFinite(d.y)&&r.push(d)}return r.length===0?null:Mt(Vt(r,2),e,i)}function Vt(n,t=0){let e=1/0,i=1/0,r=-1/0,s=-1/0;for(const o of n)e=Math.min(e,o.x),i=Math.min(i,o.y),r=Math.max(r,o.x),s=Math.max(s,o.y);return{x:e-t,y:i-t,w:r-e+t*2,h:s-i+t*2}}function Yr(n,t){let e=Math.atan2(t,n)*180/Math.PI;return e<0&&(e+=180),e}function _e(n,t){const e=n.p2.x-n.p1.x,i=n.p2.y-n.p1.y,r=Math.hypot(e,i);if(!Number.isFinite(r)||r<=1e-6)return null;const s=-i/r,o=e/r;return[{x:n.p1.x+s*t,y:n.p1.y+o*t},{x:n.p2.x+s*t,y:n.p2.y+o*t},{x:n.p2.x-s*t,y:n.p2.y-o*t},{x:n.p1.x-s*t,y:n.p1.y-o*t}]}function Vo(n,t,e,i,r,s,o,a){if(s<=0||o<=0)return null;const u=!(!!(a!=null&&a.isThreePlane)&&n.length>=t*e*3)&&((a==null?void 0:a.greenOnly)??!1),h=u?dn(a==null?void 0:a.bayerPattern,"constrained RAW edge sampling"):null,c=r.p2.x-r.p1.x,g=r.p2.y-r.p1.y,d=Math.hypot(c,g);if(!Number.isFinite(d)||d<=1e-6)return null;const f=c/d,p=g/d,m=-p,y=f,x=(r.p1.x+r.p2.x)*.5,b=(r.p1.y+r.p2.y)*.5,_=Math.abs(f)>=Math.abs(p)?1:2,k=u?Fs(n,t,e,i,0,0,h,a==null?void 0:a.greenPhase,a==null?void 0:a.blackLevel):As(n,t,e,i,{...a,globalX:0,globalY:0});if(k.length===0)return null;if(!(a!=null&&a.disableQuadraticProjection)){const P=Rs(k,a==null?void 0:a.quadraticFitPoints,r,s,o,a,!1);if(P)return P}const v=[],w=[];for(const P of k){const A=P.x-x,F=P.y-b,C=A*f+F*p;if(Math.abs(C)>s)continue;const M=A*m+F*y;Math.abs(M)>o||(v.push(M),w.push(P.value))}return v.length<8?null:pi(v,w,o*2,a==null?void 0:a.manualBinSize,_,a==null?void 0:a.preferAutoPerEdgeBin)}function Xo(n,t,e){const i=[...n],r=new Array(n.length).fill(0),s=[0,0,0];let o=-1,a=1,l=-1;for(let c=1;c<n.length-1;c++){let g=0;if(n[c]>1e-4){g=Math.atan2(e[c]*o,t[c]*o);let d=0;for(let f=-5;f<=5;f++)Math.abs(g+f*2*Math.PI-s[1])<Math.abs(g+d*2*Math.PI-s[1])&&(d=f);g+=d*2*Math.PI}c>3&&Math.abs(g-s[0])>Math.PI/2&&l<c-1&&n[c]<.5&&(a*=-1,l=c),i[c]*=a,o*=-1,s[0]=s[1],s[1]=g,s[2]=g}const u=[-.086,.343,.486,.343,-.086];for(let c=0;c<n.length-3;c++){let g=0;for(let d=-2;d<=2;d++)g+=i[Math.abs(c+d)]*u[d+2];r[c]=g}for(let c=0;c<n.length-3;c++)i[c]=r[c];const h=7;for(let c=0;c<3;c++){r.fill(0);for(let d=0;d<n.length-h;d++)if(d<h)r[d]=i[d];else{const f=Math.min(5,Math.floor((d-5)/3)),p=po[f];let m=0;for(let y=-h;y<=h;y++)m+=i[d+y]*p[y+h];r[d]=m}for(let d=n.length-h-2;d<n.length;d++)r[d]=i[d];const g=Math.abs(r[0])>1e-9?r[0]:1;for(let d=0;d<n.length;d++)i[d]=r[d]/g}for(let c=0;c<n.length;c++)i[c]=Math.abs(i[c]);return i}function Go(n,t){const e=[[0,0,0],[0,0,0],[0,0,0]],i=[0,0,0];for(let o=0;o<n.length;o++){const a=n[o],l=-t+o,u=[1,a,a*a];for(let h=0;h<3;h++){i[h]+=u[h]*l;for(let c=0;c<3;c++)e[h][c]+=u[h]*u[c]}}const r=Hs(e);if(!r)return null;const s=Qs(r,i);return[s[0],s[1],s[2]]}function zo(n,t){let e=0,i=1,r=0,s=!1,o=0;const a=Math.min(n.length,ae/16*2);for(let l=0;l<a&&!s;l++){const u=n[l];if(i>.5&&u<=.5){const h=-(u-i)*ae;Math.abs(h)>1e-9&&(r=-((.5-i-h*e)/h),o=l,s=!0)}i=u,e=l/ae}if(!s)return null;if(o>=5&&o<a-10){const l=Math.min(Math.max(2,o-9),9),u=Go(n.slice(o-l,o+l+1),l);if(u){const c=(u[0]+.5*u[1]+.25*u[2]+o)/ae;if(o>9)r=c;else{const d=(o-5)/ae/8;r=(1-d)*r+d*c}}}return r*Et*t}function Yo(n,t){if(n.length===0)return null;const e=ae,i=ae/16*4,r=new In(e),s=1,o=yo(),a=new Float32Array(501);for(let M=0;M<=500;M++)a[M]=M/500*s*2;const l=new Array(i).fill(0).map((M,S)=>S/e*s*Et),u=new Float32Array(i).fill(0),h=new Float32Array(i).fill(0);let c=0,g=[],d=[],f=[],p=[],m=[],y=[],x=[],b=null,_=0;for(const M of n){const S=M.mtfmapperOrderedDists&&M.mtfmapperOrderedVals&&M.mtfmapperOrderedDists.length===M.mtfmapperOrderedVals.length?Ss(M.mtfmapperOrderedDists,M.mtfmapperOrderedVals,M.mtfmapperEffectiveMaxDot??Lt):null,I=(S==null?void 0:S.lsfFull)??M.lsfFull,N=(S==null?void 0:S.esf)??M.esf;if(I.length<e)continue;const E=new Float32Array(e);for(let L=0;L<e;L++)E[L]=I[L]??0;r.transform(E);const X=Math.max(1e-9,Math.abs(r._real[0])),U=new Array(i).fill(0);for(let L=1;L<i;L++)U[L]=Math.atan2(r._imag[L],r._real[L]);for(let L=0;L<i;L++)u[L]+=r._real[L]/X,h[L]+=r._imag[L]/X;if(c++,_+=M.shortSidePx*.5,g.length===0){g=[...I],d=[...N];const L=new Array(i).fill(0);L[0]=1;for(let V=1;V<i;V++)L[V]=Math.hypot(r._real[V]/X,r._imag[V]/X);const B=l.map(V=>V),O=(Number.isFinite(M.zeroIndex)?M.zeroIndex:0)*(M.binSize??fi),z=vs(U,B,L,Number.POSITIVE_INFINITY,O),J=Ps(z,l,a);p=J.ptfRaw,m=J.ptfUnwrapped,y=J.ptfLinear,x=J.ptfResidual,f=J.ptfResidual,b=z.fit}}if(c===0)return null;const k=new Float32Array(i),v=new Float32Array(i),w=new Array(i).fill(0);w[0]=1;for(let M=0;M<i;M++)k[M]=u[M]/c,v[M]=h[M]/c,M>0&&(w[M]=Math.hypot(k[M],v[M]));const P=Xo(w,k,v),A=new Array(i).fill(0);for(let M=0;M<i;M++)A[M]=P[M]/o[M];const F=Array.from(a,M=>nn(M,l,A)),C=zo(A,s);return{esf:d,lsf:[],lsfCropped:g,mtf:F,ptf:f,ptfRaw:p,ptfUnwrapped:m,ptfLinear:y,ptfResidual:x,ptfPhaseFit:b,freqs:Array.from(a),mtf50:C,calcRadius:_/c}}function Wo(n,t,e,i=!1,r=0,s=!1){if(n.length===0)return null;if(n.every(w=>w.mtfmapperLike))return Yo(n);const o=4096,a=new In(o),l=1,u=new Float32Array(501);for(let w=0;w<=500;w++)u[w]=w/500*l*2;const h=new Float32Array(501).fill(0);let c=0,g=[],d=[],f=0,p=[],m=[],y=[],x=[],b=[],_=null;for(const w of n){let P=[...w.lsfFull];const A=w.binSize,F=Number.isFinite(w.zeroIndex)?w.zeroIndex:P.length/2,C=Math.max(1,Math.round((w.shortSidePx??0)*.5/Math.max(A,1e-6)));let{peakPos:M,peakIdx:S,peakVal:I}=si(P,F,C);const N=I*.2;let E=0,X=P.length-1;for(let T=S;T>=0;T--)if(P[T]<N){E=T;break}for(let T=S;T<P.length;T++)if(P[T]<N){X=T;break}const U=X-E;let L=!1;if(i&&U>0){const T=U*4,Q=[],G=[];if(r>0){const Y=Math.max(0,S-T-r),et=Math.max(0,S-T);for(let ct=Y;ct<et;ct++)Q.push(ct),G.push(P[ct]);const nt=Math.min(P.length,S+T),it=Math.min(P.length,S+T+r);for(let ct=nt;ct<it;ct++)Q.push(ct),G.push(P[ct])}else{for(let Y=0;Y<Math.max(0,S-T);Y++)Q.push(Y),G.push(P[Y]);for(let Y=Math.min(P.length,S+T);Y<P.length;Y++)Q.push(Y),G.push(P[Y])}if(Q.length>2){const{slope:Y,intercept:et}=ir(Q,G);for(let nt=0;nt<P.length;nt++)P[nt]=P[nt]-(Y*nt+et);({peakPos:M,peakIdx:S,peakVal:I}=si(P,F,C)),L=!0}}let B=0,O=0;if(t>0)O=t,B=Math.round(t/A);else{const T=I*.2;let Q=0,G=P.length-1;for(let it=S;it>=0;it--)if(P[it]<T){Q=it;break}for(let it=S;it<P.length;it++)if(P[it]<T){G=it;break}const et=(G-Q)*A;let nt=Math.max(2,et*8);O=nt,B=Math.round(nt/A)}f+=O;const z=Math.max(0,Math.floor(F-B)),J=Math.min(P.length,Math.ceil(F+B)),V=P.slice(z,J);if(V.length<8)continue;const q=new Float32Array(o).fill(0),W=new Array(V.length).fill(0);for(let T=0;T<V.length;T++){let Q=1;s&&(Q=.5*(1-Math.cos(2*Math.PI*T/(V.length-1)))),W[T]=V[T]*Q}const tt=Math.max(0,Math.min(V.length-1,M-z));for(let T=0;T<o;T++)q[T]=ko(W,o,T+tt);a.transform(q);const at=[],D=[],Z=[];for(let T=0;T<=o/2;T++){const Q=a._real[T],G=a._imag[T],Y=Math.sqrt(Q*Q+G*G);at.push(Y),D.push(T/(o*A)*l),Z.push(Math.atan2(G,Q))}const H=at[0];if(H>0){for(let T=0;T<=500;T++){const Q=u[T],Y=go(Q,A);h[T]+=nn(Q,D,at)/H/Y}if(c++,g.length===0){g=Fo(V,tt,(V.length-1)/2),d=L?jo(P):w.esf;const T=at.map(et=>et/H),Q=D.map(et=>et),G=vs(Z,Q,T,Number.POSITIVE_INFINITY,0),Y=Ps(G,D,u);m=Y.ptfRaw,y=Y.ptfUnwrapped,x=Y.ptfLinear,b=Y.ptfResidual,p=Y.ptfResidual,_=G.fit}}}if(c===0)return null;const k=Array.from(h).map(w=>w/c);let v=null;for(let w=0;w<k.length-1;w++)if(k[w]>=.5&&k[w+1]<.5){v=u[w]+(.5-k[w])*(u[w+1]-u[w])/(k[w+1]-k[w]);break}return{esf:d,lsf:[],lsfCropped:g,mtf:k,ptf:p,ptfRaw:m,ptfUnwrapped:y,ptfLinear:x,ptfResidual:b,ptfPhaseFit:_,freqs:Array.from(u),mtf50:v,calcRadius:f/c}}function jo(n){const t=new Array(n.length).fill(0);let e=0;for(let i=0;i<n.length;i++)e+=n[i],t[i]=e;return t}function nn(n,t,e){if(n<=t[0])return e[0];if(n>=t[t.length-1])return e[e.length-1];let i=0;for(;n>t[i+1];)i++;const r=(n-t[i])/(t[i+1]-t[i]);return e[i]+r*(e[i+1]-e[i])}function ar(n){return{...ji,...n,gradientPercentiles:n!=null&&n.gradientPercentiles&&n.gradientPercentiles.length>0?n.gradientPercentiles:ji.gradientPercentiles}}function Ho(n){return!n||n.length===0?void 0:[Number.isFinite(n[0])?n[0]:0,Number.isFinite(n[1])?n[1]:Number.isFinite(n[0])?n[0]:0,Number.isFinite(n[2])?n[2]:Number.isFinite(n[0])?n[0]:0,Number.isFinite(n[3])?n[3]:Number.isFinite(n[0])?n[0]:0]}function Qo(n,t){const e=n.width,i=n.height,r=n.data,s=dn(n.bayerPattern,"RAW SFR detection"),o=Ho(n.blackLevels),a=new Float32Array(e*i),l=(_,k)=>_<0||k<0||_>=e||k>=i?null:Math.max(0,r[k*e+_]-Mn(o,_,k));let u=1/0,h=-1/0;for(let _=0;_<i;_++){const k=_*e;for(let v=0;v<e;v++){const w=k+v;let P=0;if(vt(v,_,s,t))P=l(v,_)??0;else{const A=[],F=l(v-1,_),C=l(v+1,_),M=l(v,_-1),S=l(v,_+1);if(F!==null&&vt(v-1,_,s,t)&&A.push(F),C!==null&&vt(v+1,_,s,t)&&A.push(C),M!==null&&vt(v,_-1,s,t)&&A.push(M),S!==null&&vt(v,_+1,s,t)&&A.push(S),A.length>0)P=ri(A);else{const I=[],N=l(v-1,_-1),E=l(v+1,_-1),X=l(v-1,_+1),U=l(v+1,_+1);N!==null&&vt(v-1,_-1,s,t)&&I.push(N),E!==null&&vt(v+1,_-1,s,t)&&I.push(E),X!==null&&vt(v-1,_+1,s,t)&&I.push(X),U!==null&&vt(v+1,_+1,s,t)&&I.push(U),P=ri(I)}}a[w]=P,P<u&&(u=P),P>h&&(h=P)}}if(!Number.isFinite(u)||!Number.isFinite(h)||h<=u+1e-9)return new Uint8Array(e*i);const c=1024,g=new Uint32Array(c),d=h-u;for(let _=0;_<a.length;_++){const k=Math.max(0,Math.min(1,(a[_]-u)/d)),v=Math.min(c-1,Math.max(0,Math.floor(k*(c-1))));g[v]++}const f=a.length,p=_=>{const k=f*_;let v=0;for(let w=0;w<c;w++)if(v+=g[w],v>=k)return u+w/Math.max(1,c-1)*d;return h},m=p(.01),y=p(.99),x=Math.max(1e-9,y-m),b=new Uint8Array(e*i);for(let _=0;_<a.length;_++){const k=Math.max(0,Math.min(1,(a[_]-m)/x));b[_]=Math.round(k*255)}return b}function qo(n,t,e){const i=new Float32Array(n.length),r=new Float32Array(n.length),s=new Float32Array(n.length);for(let o=1;o<e-1;o++)for(let a=1;a<t-1;a++){const l=o*t+a,u=n[(o-1)*t+(a-1)],h=n[(o-1)*t+a],c=n[(o-1)*t+(a+1)],g=n[o*t+(a-1)],d=n[o*t+(a+1)],f=n[(o+1)*t+(a-1)],p=n[(o+1)*t+a],m=n[(o+1)*t+(a+1)],y=-u-2*g-f+(c+2*d+m),x=-u-2*h-c+(f+2*p+m);i[l]=y,r[l]=x,s[l]=Math.hypot(y,x)}return{gx:i,gy:r,magnitude:s}}function Ko(n,t){let e=0,i=0;for(let l=0;l<n.length;l++){const u=n[l];!Number.isFinite(u)||u<=1e-6||(e=Math.max(e,u),i++)}if(i===0||e<=1e-6)return[];const r=1024,s=new Uint32Array(r);for(let l=0;l<n.length;l++){const u=n[l];if(!Number.isFinite(u)||u<=1e-6)continue;const h=Math.max(0,Math.min(1,u/e)),c=Math.min(r-1,Math.floor(h*(r-1)));s[c]++}const o=t&&t.length>0?t:ji.gradientPercentiles,a=[];for(const l of o){const u=i*l;let h=0;for(let c=0;c<r;c++)if(h+=s[c],h>=u){a.push(c/Math.max(1,r-1)*e);break}}return Array.from(new Set(a.filter(l=>l>0))).sort((l,u)=>u-l)}function $o(n,t){const e=new Uint8Array(n.length);for(let i=0;i<n.length;i++)e[i]=n[i]>=t?1:0;return e}const Jo=256*256;function Zo(n,t,e){if(n.length>=Jo){const s=uo.compute(n,t,e);if(s)return{gray:s.blurredGray,gradient:{gx:s.gx,gy:s.gy,magnitude:s.magnitude},backend:"webgl"}}const i=cl(n,t,e),r=qo(i,t,e);return{gray:i,gradient:r,backend:"cpu"}}function tl(n,t,e,i){let r=n;for(let s=0;s<i;s++){const o=new Uint8Array(n.length);for(let a=0;a<e;a++)for(let l=0;l<t;l++){let u=0;for(let h=-1;h<=1&&!u;h++){const c=a+h;if(!(c<0||c>=e))for(let g=-1;g<=1;g++){const d=l+g;if(!(d<0||d>=t)&&r[c*t+d]){u=1;break}}}o[a*t+l]=u}r=o}return r}function el(n,t,e){const i=new Int32Array(n.length),r=[];let s=1;for(let o=0;o<n.length;o++){if(!n[o]||i[o]!==0)continue;const a=[o];i[o]=s;let l=0,u=t,h=e,c=0,g=0,d=0,f=!1;for(;l<a.length;){const p=a[l++],m=p%t,y=Math.floor(p/t);d++,u=Math.min(u,m),h=Math.min(h,y),c=Math.max(c,m),g=Math.max(g,y),(m===0||y===0||m===t-1||y===e-1)&&(f=!0);for(let x=-1;x<=1;x++)for(let b=-1;b<=1;b++){if(b===0&&x===0)continue;const _=m+b,k=y+x;if(_<0||k<0||_>=t||k>=e)continue;const v=k*t+_;!n[v]||i[v]!==0||(i[v]=s,a.push(v))}}r.push({label:s,x:u,y:h,w:c-u+1,h:g-h+1,area:d,touchesBorder:f}),s++}return{labels:i,components:r}}function Ns(n,t){const e=Math.hypot(n,t);if(!Number.isFinite(e)||e<=1e-9)return null;let i=n/e,r=t/e;return(i<0||Math.abs(i)<=1e-9&&r<0)&&(i=-i,r=-r),{x:i,y:r}}function Ie(n,t){if(n.length===0)return 0;const e=[...n].sort((o,a)=>o.value-a.value),i=e.reduce((o,a)=>o+Math.max(0,a.weight),0);if(i<=0)return e[Math.floor((e.length-1)*t)].value;const r=Math.max(0,Math.min(1,t))*i;let s=0;for(const o of e)if(s+=Math.max(0,o.weight),s>=r)return o.value;return e[e.length-1].value}function Wr(n){const t=n.filter(i=>Number.isFinite(i)).sort((i,r)=>i-r);if(t.length===0)return 0;const e=i=>{if(i.length===1)return i[0];if(i.length===2)return(i[0]+i[1])*.5;const r=Math.ceil(i.length*.5);let s=0,o=1/0;for(let a=0;a+r-1<i.length;a++){const l=i[a+r-1]-i[a];l<o&&(o=l,s=a)}return e(i.slice(s,s+r))};return e(t)}function nl(n,t,e,i,r,s,o){const a=[];for(let l=r.y;l<r.y+r.h;l++)for(let u=r.x;u<r.x+r.w;u++){const h=l*s+u;if(n[h]!==t||!e[h])continue;const c=i.magnitude[h];!Number.isFinite(c)||c<=1e-6||a.push({x:u,y:l,weight:c,gx:i.gx[h],gy:i.gy[h]})}return a}function il(n){let t=0,e=0,i=0,r=0,s=0;for(const l of n){t+=l.weight,e+=l.x*l.weight,i+=l.y*l.weight;const u=Math.hypot(l.gx,l.gy);if(!Number.isFinite(u)||u<=1e-6)continue;const h=-l.gy/u,c=l.gx/u;r+=l.weight*(h*h-c*c),s+=l.weight*(2*h*c)}if(t<=0)return null;e/=t,i/=t;const o=.5*Math.atan2(s,r),a=Ns(Math.cos(o),Math.sin(o));return a?{centerX:e,centerY:i,dirX:a.x,dirY:a.y,orthoX:-a.y,orthoY:a.x}:null}function Qn(n,t){let e=0,i=0;const r=-t.dirY,s=t.dirX;for(const o of n){const a=(o.x-t.pointX)*r+(o.y-t.pointY)*s;i+=o.weight*a*a,e+=o.weight}return e<=0?1/0:Math.sqrt(i/e)}function jr(n,t,e,i,r){const s=Math.max(0,Math.min(t-1,i)),o=Math.max(0,Math.min(e-1,r)),a=Math.floor(s),l=Math.floor(o),u=Math.min(t-1,a+1),h=Math.min(e-1,l+1),c=s-a,g=o-l,d=n[l*t+a],f=n[l*t+u],p=n[h*t+a],m=n[h*t+u],y=d+(f-d)*c,x=p+(m-p)*c;return y+(x-y)*g}function rl(n,t,e,i,r,s,o,a){const l=jr(n,t,e,i-s*a,r-o*a);return jr(n,t,e,i+s*a,r+o*a)-l}function qn(n,t,e,i,r){const s=Math.max(1e-6,e-t);if(n.length===0||!Number.isFinite(s))return{points:[],coverageRatio:0,centerCoverageRatio:0};const o=Math.max(1.5,Math.min(4,s/18)),a=Math.max(1,Math.ceil(s/o)),l=new Map;for(const f of n){const p=i(f);if(!Number.isFinite(p)||p<t||p>e)continue;const m=Math.max(0,Math.min(a-1,Math.floor((p-t)/o))),y=f.weight/(1+Math.abs(r(f))),x=l.get(m);(!x||y>x.score)&&l.set(m,{point:f,score:y})}const u=Array.from(l.values()).sort((f,p)=>i(f.point)-i(p.point)).map(f=>f.point),h=Math.max(0,Math.floor(a*.3)),c=Math.max(h+1,Math.ceil(a*.7));let g=0;for(let f=h;f<c;f++)l.has(f)&&g++;const d=Math.max(1,c-h);return{points:u,coverageRatio:u.length/a,centerCoverageRatio:g/d}}function Kn(n,t){const e=n.dirX*t.dirY-n.dirY*t.dirX;if(!Number.isFinite(e)||Math.abs(e)<=1e-6)return null;const i=t.pointX-n.pointX,r=t.pointY-n.pointY,s=(i*t.dirY-r*t.dirX)/e;return{x:n.pointX+n.dirX*s,y:n.pointY+n.dirY*s}}function sl(n){if(n.length<3)return 0;let t=0;for(let e=0;e<n.length;e++){const i=n[e],r=n[(e+1)%n.length];t+=i.x*r.y-r.x*i.y}return t*.5}function al(n,t,e,i,r,s,o,a,l){const u=nl(i,r.label,s,o,r,t),h=u.map(R=>({x:R.x,y:R.y}));if(u.length<l.minEdgePoints)return{candidate:null,failureStage:"min_edge_points",pointsCount:u.length,strongEdgePoints:h};const c=il(u);if(!c)return{candidate:null,failureStage:"dominant_axes",pointsCount:u.length,strongEdgePoints:h};const g=u.map(R=>{const lt=R.x-c.centerX,ht=R.y-c.centerY;return{...R,u:lt*c.dirX+ht*c.dirY,v:lt*c.orthoX+ht*c.orthoY}}),d={x:c.centerX,y:c.centerY},f=Ie(g.map(R=>({value:R.u,weight:R.weight})),l.extentQuantileLow),p=Ie(g.map(R=>({value:R.u,weight:R.weight})),l.extentQuantileHigh),m=Ie(g.map(R=>({value:R.v,weight:R.weight})),l.extentQuantileLow),y=Ie(g.map(R=>({value:R.v,weight:R.weight})),l.extentQuantileHigh),x=Math.max(1e-6,Math.max(Math.abs(f),Math.abs(p))),b=Math.max(1e-6,Math.max(Math.abs(m),Math.abs(y))),_=72,k=360/_,v=Array.from({length:_},()=>[]),w=R=>{let lt=R%360;return lt<0&&(lt+=360),lt},P=(R,lt)=>{const ht=Math.abs(w(R)-w(lt));return Math.min(ht,360-ht)};g.forEach(R=>{const lt=R.u/x,ht=R.v/b,At=w(Math.atan2(ht,lt)*180/Math.PI),Yt=Math.hypot(lt,ht),Nt=Math.max(0,Math.min(_-1,Math.floor(At/k)));v[Nt].push({point:R,angleDeg:At,normRadius:Yt})});const A=v.map(R=>R.length>0?Wr(R.map(lt=>lt.normRadius)):-1/0),F=(R,lt)=>{let ht=-1,At=-1/0;for(let yt=0;yt<v.length;yt++){if(v[yt].length===0)continue;const me=(yt+.5)*k;if(P(me,R)>45||lt.some(Ys=>P(me,Ys)<45))continue;const Te=A[yt];Te>At&&(At=Te,ht=yt)}let Yt=ht>=0?(ht+.5)*k:R,Nt=ht>=0?v[ht]:g.map(yt=>{const pe=yt.u/x,me=yt.v/b;return{point:yt,angleDeg:w(Math.atan2(me,pe)*180/Math.PI),normRadius:Math.hypot(pe,me)}}).filter(yt=>P(yt.angleDeg,R)<=45&&!lt.some(pe=>P(yt.angleDeg,pe)<45));if(Nt.length===0&&(Nt=g.map(yt=>{const pe=yt.u/x,me=yt.v/b;return{point:yt,angleDeg:w(Math.atan2(me,pe)*180/Math.PI),normRadius:Math.hypot(pe,me)}}).filter(yt=>P(yt.angleDeg,R)<=45),Yt=R),Nt.length===0)return{x:g[0].x,y:g[0].y,u:g[0].u,v:g[0].v,angleDeg:R};const zn=ht>=0?A[ht]:Wr(Nt.map(yt=>yt.normRadius));let Se=0,Pn=0,mr=0,gr=0,yr=0;for(const yt of Nt){const pe=P(yt.angleDeg,R)/45,me=Math.abs(yt.normRadius-zn),Te=Math.max(1e-6,yt.point.weight)/(1+pe*2+me*6);Se+=Te,Pn+=yt.point.x*Te,mr+=yt.point.y*Te,gr+=yt.point.u*Te,yr+=yt.point.v*Te}return Se>0?{x:Pn/Se,y:mr/Se,u:gr/Se,v:yr/Se,angleDeg:Yt}:{x:Nt[0].point.x,y:Nt[0].point.y,u:Nt[0].point.u,v:Nt[0].point.v,angleDeg:Nt[0].angleDeg}},C=F(225,[]),M=F(315,[C.angleDeg]),S=F(45,[C.angleDeg,M.angleDeg]),I=F(135,[C.angleDeg,M.angleDeg,S.angleDeg]),N=[{x:C.x,y:C.y},{x:M.x,y:M.y},{x:S.x,y:S.y},{x:I.x,y:I.y}],E=p-f,X=y-m,U=Math.min(E,X),L=Math.max(E,X);if(!Number.isFinite(U)||U<l.minSpanPx||L/Math.max(1,U)>l.maxAspectRatio)return{candidate:null,failureStage:"span_aspect",pointsCount:u.length,minSpan:U,maxSpan:L,axisCentroid:d,axisExtremePoints:N,strongEdgePoints:h};const B=Math.max(l.bandMinPx,Math.min(l.bandMaxPx,U*l.bandScale)),O=Math.max(1,Math.min(3,B*.55)),z=Math.max(a,0),J=void 0,V=void 0,q=R=>R.map(lt=>({x:lt.x,y:lt.y,weight:lt.weight})),W=R=>R.map(lt=>({x:lt.x,y:lt.y})),tt=(R,lt,ht)=>R.filter(At=>{if(!Number.isFinite(At.weight)||At.weight<z)return!1;const Yt=rl(n,t,e,At.x,At.y,lt,ht,O);return Number.isFinite(Yt)&&Yt>=l.minPointContrast}),at=f,D=p,Z=m,H=y,T=l.minCoverageRatio,Q=l.minCenterCoverageRatio,G=[],Y=[],et=[],nt=[],it=[],ct=(R,lt,ht,At,Yt,Nt)=>(ht-R)*(Nt-lt)-(At-lt)*(Yt-R),dt=R=>R>1e-6?1:R<-1e-6?-1:0,rt=[{u:(C.u+M.u)*.5,v:(C.v+M.v)*.5},{u:(M.u+S.u)*.5,v:(M.v+S.v)*.5},{u:(S.u+I.u)*.5,v:(S.v+I.v)*.5},{u:(I.u+C.u)*.5,v:(I.v+C.v)*.5}],Bt=(R,lt)=>{const ht=dt(ct(C.u,C.v,S.u,S.v,R,lt)),At=dt(ct(M.u,M.v,I.u,I.v,R,lt));return`${ht},${At}`},Xt=new Map;rt.forEach((R,lt)=>{Xt.set(Bt(R.u,R.v),lt)});for(const R of g){if(!Number.isFinite(R.u)||!Number.isFinite(R.v)){it.push(R);continue}let ht=Xt.get(Bt(R.u,R.v))??-1;if(ht<0){let At=1/0;for(let Yt=0;Yt<rt.length;Yt++){const Nt=rt[Yt],zn=(R.u-Nt.u)/x,Se=(R.v-Nt.v)/b,Pn=zn*zn+Se*Se;Pn<At&&(At=Pn,ht=Yt)}}ht===0?G.push(R):ht===1?Y.push(R):ht===2?et.push(R):ht===3?nt.push(R):it.push(R)}const j=[...G,...et],ut=[...Y,...nt],ot={dir:j.length,ortho:ut.length,unassigned:g.length-j.length-ut.length},pt=G.length>=l.minSidePoints?Ie(G.map(R=>({value:R.v,weight:R.weight})),.5):m,mt=et.length>=l.minSidePoints?Ie(et.map(R=>({value:R.v,weight:R.weight})),.5):y,ee=nt.length>=l.minSidePoints?Ie(nt.map(R=>({value:R.u,weight:R.weight})),.5):f,Dt=Y.length>=l.minSidePoints?Ie(Y.map(R=>({value:R.u,weight:R.weight})),.5):p,_t=[{x:(C.x+M.x)*.5,y:(C.y+M.y)*.5},{x:(M.x+S.x)*.5,y:(M.y+S.y)*.5},{x:(S.x+I.x)*.5,y:(S.y+I.y)*.5},{x:(I.x+C.x)*.5,y:(I.y+C.y)*.5}],ce=G.filter(R=>Math.abs(R.v-pt)<=B&&R.u>=at&&R.u<=D),$t=et.filter(R=>Math.abs(R.v-mt)<=B&&R.u>=at&&R.u<=D),Ae=nt.filter(R=>Math.abs(R.u-ee)<=B&&R.v>=Z&&R.v<=H),$=Y.filter(R=>Math.abs(R.u-Dt)<=B&&R.v>=Z&&R.v<=H),Ft=[ce.length,$.length,$t.length,Ae.length],Tt=[W(ce),W($),W($t),W(Ae)],xt=tt(ce,-c.orthoX,-c.orthoY),wt=tt($t,c.orthoX,c.orthoY),ne=tt(Ae,-c.dirX,-c.dirY),Jt=tt($,c.dirX,c.dirY),Sn=[xt.length,Jt.length,wt.length,ne.length],Oe=[W(xt),W(Jt),W(wt),W(ne)],Ve=qn(xt,at,D,R=>R.u,R=>R.v-pt),Xe=qn(Jt,Z,H,R=>R.v,R=>R.u-Dt),Ge=qn(wt,at,D,R=>R.u,R=>R.v-mt),ze=qn(ne,Z,H,R=>R.v,R=>R.u-ee),Nn=(R,lt)=>R.slice().sort((ht,At)=>lt(ht)-lt(At)),Ln=Nn(xt,R=>R.u),En=Nn(Jt,R=>R.v),Un=Nn(wt,R=>R.u),Bn=Nn(ne,R=>R.v),mi=[Ln.length,En.length,Un.length,Bn.length],Ht=[Ve.coverageRatio,Xe.coverageRatio,Ge.coverageRatio,ze.coverageRatio];Ve.centerCoverageRatio,Xe.centerCoverageRatio,Ge.centerCoverageRatio,ze.centerCoverageRatio;const Gt=[W(Ln),W(En),W(Un),W(Bn)],Zt={axisPointCounts:ot,sideBandPointCounts:Ft,sideContrastPointCounts:Sn,gradientThreshold:a,pointAxisMinDot:J,pointAxisMargin:V,bandWidth:B,minPointContrast:l.minPointContrast,minCoverageRatio:T,minCenterCoverageRatio:Q,axisCentroid:d,axisExtremePoints:N,axisSideCenters:_t,strongEdgePoints:h,axisDirPoints:W(j),axisOrthoPoints:W(ut),axisUnassignedPoints:W(it),sideBandPoints:Tt,sideContrastPoints:Oe};if(Ln.length<l.minSidePoints||Un.length<l.minSidePoints||Bn.length<l.minSidePoints||En.length<l.minSidePoints)return{candidate:null,failureStage:"min_side_points",pointsCount:u.length,minSpan:U,maxSpan:L,sidePointCounts:mi,sideCoverageRatios:Ht,...Zt,sideFitPoints:Gt};if(Ve.coverageRatio<T||Xe.coverageRatio<T||Ge.coverageRatio<T||ze.coverageRatio<T||Ve.centerCoverageRatio<Q||Xe.centerCoverageRatio<Q||Ge.centerCoverageRatio<Q||ze.centerCoverageRatio<Q)return{candidate:null,failureStage:"side_coverage",pointsCount:u.length,minSpan:U,maxSpan:L,sidePointCounts:mi,sideCoverageRatios:Ht,...Zt,sideFitPoints:Gt};const gi=Ln,yi=Un,xi=Bn,bi=En,ie=mi,Ye=Ce(q(gi)),We=Ce(q(yi)),je=Ce(q(xi)),He=Ce(q(bi));if(!Ye||!We||!je||!He)return{candidate:null,failureStage:"fit_lines",pointsCount:u.length,minSpan:U,maxSpan:L,sidePointCounts:ie,sideCoverageRatios:Ht,...Zt,sideFitPoints:Gt};const we=Uo(q(gi),Ye,q(bi),He,q(yi),We,q(xi),je),Dn=l.minAxisDot,On=(R,lt,ht)=>Math.abs(R.dirX*lt+R.dirY*ht),zt=[On(Ye,c.dirX,c.dirY),On(He,c.orthoX,c.orthoY),On(We,c.dirX,c.dirY),On(je,c.orthoX,c.orthoY)];if(zt[0]<Dn||zt[1]<Dn||zt[2]<Dn||zt[3]<Dn)return{candidate:null,failureStage:"axis_alignment",pointsCount:u.length,minSpan:U,maxSpan:L,sidePointCounts:ie,sideCoverageRatios:Ht,axisDots:zt,...Zt,sideFitPoints:Gt,sideFitLines:we};const re=Math.max(l.residualLimitFloor,B*l.residualLimitScale),gt=[Qn(q(gi),Ye),Qn(q(yi),We),Qn(q(xi),je),Qn(q(bi),He)],_i=[gt[0],gt[3],gt[1],gt[2]],Mi=Math.max(...gt),ue=Kn(Ye,je),he=Kn(Ye,He),de=Kn(We,He),fe=Kn(We,je);if(!ue||!he||!de||!fe)return{candidate:null,failureStage:"corners",pointsCount:u.length,minSpan:U,maxSpan:L,sidePointCounts:ie,sideCoverageRatios:Ht,axisDots:zt,sideResiduals:[gt[0],gt[3],gt[1],gt[2]],residualLimit:re,...Zt,sideFitPoints:Gt,sideFitLines:we};const vn=[ue,he,de,fe],se=Math.abs(sl(vn));if(!Number.isFinite(se)||se<l.minQuadArea)return{candidate:null,failureStage:"quad_area",pointsCount:u.length,minSpan:U,maxSpan:L,sidePointCounts:ie,sideCoverageRatios:Ht,axisDots:zt,sideResiduals:[gt[0],gt[3],gt[1],gt[2]],residualLimit:re,...Zt,sideFitPoints:Gt,sideFitLines:we,quadArea:se};const wi=Math.hypot(he.x-ue.x,he.y-ue.y),Si=Math.hypot(de.x-he.x,de.y-he.y),vi=Math.hypot(de.x-fe.x,de.y-fe.y),Pi=Math.hypot(fe.x-ue.x,fe.y-ue.y),Qe=[wi,Si,vi,Pi],Ci=Math.min(wi,Si,vi,Pi),Vs=Math.max(wi,Si,vi,Pi);if(!Number.isFinite(Ci)||Ci<l.minSideLength||Vs/Math.max(1,Ci)>l.maxAspectRatio)return{candidate:null,failureStage:"side_length",pointsCount:u.length,minSpan:U,maxSpan:L,sidePointCounts:ie,sideCoverageRatios:Ht,axisDots:zt,sideResiduals:[gt[0],gt[3],gt[1],gt[2]],residualLimit:re,...Zt,sideFitPoints:Gt,sideFitLines:we,quadArea:se,sideLengths:Qe};const qe=Ns(he.x-ue.x+(de.x-fe.x),he.y-ue.y+(de.y-fe.y));if(!qe)return{candidate:null,failureStage:"corners",pointsCount:u.length,minSpan:U,maxSpan:L,sidePointCounts:ie,sideCoverageRatios:Ht,axisDots:zt,sideResiduals:[gt[0],gt[3],gt[1],gt[2]],residualLimit:re,...Zt,sideFitPoints:Gt,sideFitLines:we,quadArea:se,sideLengths:Qe};const dr={x:-qe.y,y:qe.x},ki=(ue.x+he.x+de.x+fe.x)*.25,Fi=(ue.y+he.y+de.y+fe.y)*.25,Vn=vn.map(R=>{const lt=R.x-ki,ht=R.y-Fi;return{u:lt*qe.x+ht*qe.y,v:lt*dr.x+ht*dr.y}}),Xn=(Math.max(...Vn.map(R=>R.u))-Math.min(...Vn.map(R=>R.u)))*.5,Gn=(Math.max(...Vn.map(R=>R.v))-Math.min(...Vn.map(R=>R.v)))*.5;if(!Number.isFinite(Xn)||!Number.isFinite(Gn)||Math.min(Xn,Gn)<6)return{candidate:null,failureStage:"box_size",pointsCount:u.length,minSpan:U,maxSpan:L,sidePointCounts:ie,sideCoverageRatios:Ht,axisDots:zt,sideResiduals:[gt[0],gt[3],gt[1],gt[2]],residualLimit:re,sideFitPoints:Gt,quadArea:se,sideLengths:Qe};const K=fl(n,t,e,vn,Ye,He,We,je,ki,Fi,Xn,Gn,l.innerPurityStdScale,l.outerMeanSpreadLimit);if(!Number.isFinite(Mi)||Mi>re)return{candidate:null,failureStage:"residual",pointsCount:u.length,minSpan:U,maxSpan:L,sidePointCounts:ie,sideCoverageRatios:Ht,axisDots:zt,sideResiduals:_i,residualLimit:re,...Zt,sideFitPoints:Gt,sideFitLines:we,quadArea:se,sideLengths:Qe,outerContrast:K.contrast,outerUniformityOk:K.ok,outerMeanSpread:K.meanSpread,outerMeanSpreadLimit:K.meanSpreadLimit,outerAvgStd:K.avgStd,outerAvgStdLimit:K.avgStdLimit,outerSideMeans:K.outerSideMeans,outerSideStds:K.outerSideStds,outerSideStdLimit:K.outerSideStdLimit,outerSideQuads:K.outerSideQuads,innerSideUniformityOk:K.innerSideOk,innerSideStds:K.innerSideStds,innerSideStdLimit:K.innerSideStdLimit,innerSideQuads:K.innerSideQuads};const fr=l.filterBlockPurity&&(!K.ok||!K.innerSideOk);if(fr||K.contrast<l.minOuterContrast)return{candidate:null,failureStage:fr?K.ok?"inner_roi_uniformity":"outer_uniformity":"outer_contrast",pointsCount:u.length,minSpan:U,maxSpan:L,sidePointCounts:ie,sideCoverageRatios:Ht,axisDots:zt,sideResiduals:_i,residualLimit:re,...Zt,sideFitPoints:Gt,sideFitLines:we,quadArea:se,sideLengths:Qe,outerContrast:K.contrast,outerUniformityOk:K.ok,outerMeanSpread:K.meanSpread,outerMeanSpreadLimit:K.meanSpreadLimit,outerAvgStd:K.avgStd,outerAvgStdLimit:K.avgStdLimit,outerSideMeans:K.outerSideMeans,outerSideStds:K.outerSideStds,outerSideStdLimit:K.outerSideStdLimit,outerSideQuads:K.outerSideQuads,innerSideUniformityOk:K.innerSideOk,innerSideStds:K.innerSideStds,innerSideStdLimit:K.innerSideStdLimit,innerSideQuads:K.innerSideQuads};const pr=Mt(Vt(vn,1),t,e);if(!pr)return{candidate:null,failureStage:"bbox",pointsCount:u.length,minSpan:U,maxSpan:L,sidePointCounts:ie,sideCoverageRatios:Ht,axisDots:zt,sideResiduals:[gt[0],gt[3],gt[1],gt[2]],residualLimit:re,...Zt,sideFitPoints:Gt,sideFitLines:we,quadArea:se,sideLengths:Qe,outerContrast:K.contrast,outerUniformityOk:K.ok,outerMeanSpread:K.meanSpread,outerMeanSpreadLimit:K.meanSpreadLimit,outerAvgStd:K.avgStd,outerAvgStdLimit:K.avgStdLimit,outerSideMeans:K.outerSideMeans,outerSideStds:K.outerSideStds,outerSideStdLimit:K.outerSideStdLimit,outerSideQuads:K.outerSideQuads,innerSideUniformityOk:K.innerSideOk,innerSideStds:K.innerSideStds,innerSideStdLimit:K.innerSideStdLimit,innerSideQuads:K.innerSideQuads};const Xs=1/(1+Mi/Math.max(1,re)),Gs=l.filterBlockPurity?K.score:1,zs=K.contrast*Gs*Xs*Math.sqrt(se);return{candidate:{centerX:ki,centerY:Fi,dirX:qe.x,dirY:qe.y,halfWidth:Xn,halfHeight:Gn,score:zs,bbox:pr,corners:vn,sideFitPoints:Gt,outerSideMeans:K.outerSideMeans,outerSideQuads:K.outerSideQuads},failureStage:null,pointsCount:u.length,minSpan:U,maxSpan:L,sidePointCounts:ie,sideCoverageRatios:Ht,axisDots:zt,sideResiduals:_i,residualLimit:re,...Zt,sideFitPoints:Gt,sideFitLines:we,quadArea:se,sideLengths:Qe,outerContrast:K.contrast,outerUniformityOk:K.ok,outerMeanSpread:K.meanSpread,outerMeanSpreadLimit:K.meanSpreadLimit,outerAvgStd:K.avgStd,outerAvgStdLimit:K.avgStdLimit,outerSideMeans:K.outerSideMeans,outerSideStds:K.outerSideStds,outerSideStdLimit:K.outerSideStdLimit,outerSideQuads:K.outerSideQuads,innerSideUniformityOk:K.innerSideOk,innerSideStds:K.innerSideStds,innerSideStdLimit:K.innerSideStdLimit,innerSideQuads:K.innerSideQuads}}function ol(n,t,e,i,r,s,o,a,l){return al(n,t,e,i,r,s,o,a,l).candidate}function ll(n,t,e,i,r,s){const o=ar(r),a=Math.max(i*8,i+128),l=ul(n,t,e,o.downsampleMaxSide);s==null||s("Detecting candidates: downsampling...",.02),s==null||s("Detecting candidates: edge stage...",.06);const u=Zo(l.gray,l.width,l.height),h=u.gray,c=u.gradient;s==null||s(`Detecting candidates: gradient (${u.backend==="webgl"?"WebGL1":"CPU"})...`,.1);const g=Ko(c.magnitude,o.gradientPercentiles),d=l.width*l.height,f=Math.max(o.minComponentAreaPx,Math.round(d*o.minComponentAreaRatio)),p=Math.max(f+1,Math.round(d*o.maxComponentAreaRatio)),m=[],y=Math.max(1,g.reduce((C,M,S)=>C+(S<=1,2),0));let x=0;for(let C=0;C<g.length;C++){const M=g[C],S=$o(c.magnitude,M),I=C<=1?[3,2]:[2,1];for(const N of I){const E=x/y;s==null||s(`Detecting candidates: threshold ${C+1}/${g.length}, dilate ${N}`,.12+.78*E);const X=tl(S,l.width,l.height,N),{labels:U,components:L}=el(X,l.width,l.height);for(const B of L){if(B.touchesBorder||B.area<f||B.area>p)continue;const O=ol(h,l.width,l.height,U,B,S,c,M,o);if(!O)continue;const z=1/l.scale,J=O.corners.map(V=>({x:V.x*z,y:V.y*z}));m.push({centerX:O.centerX*z,centerY:O.centerY*z,dirX:O.dirX,dirY:O.dirY,halfWidth:O.halfWidth*z,halfHeight:O.halfHeight*z,score:O.score,bbox:{x:O.bbox.x*z,y:O.bbox.y*z,w:O.bbox.w*z,h:O.bbox.h*z},corners:J,sideFitPoints:O.sideFitPoints?[O.sideFitPoints[0].map(V=>({x:V.x*z,y:V.y*z})),O.sideFitPoints[1].map(V=>({x:V.x*z,y:V.y*z})),O.sideFitPoints[2].map(V=>({x:V.x*z,y:V.y*z})),O.sideFitPoints[3].map(V=>({x:V.x*z,y:V.y*z}))]:void 0,outerSideMeans:O.outerSideMeans,outerSideQuads:O.outerSideQuads?[O.outerSideQuads[0].map(V=>({x:V.x*z,y:V.y*z})),O.outerSideQuads[1].map(V=>({x:V.x*z,y:V.y*z})),O.outerSideQuads[2].map(V=>({x:V.x*z,y:V.y*z})),O.outerSideQuads[3].map(V=>({x:V.x*z,y:V.y*z}))]:void 0})}m.length>a&&(m.sort((B,O)=>O.score-B.score),m.length=a),x++}}console.log(`[SFR Auto Detect] Candidate pool before dedupe: ${m.length}`),s==null||s(`Detecting candidates: deduplicating (0/${Math.max(1,Math.min(m.length,Math.max(i*4,i+32)))})...`,.94),m.sort((C,M)=>M.score-C.score);const b=Math.max(i*4,i+32),_=m.length>b?m.slice(0,b):m,k=[];if(_.length<=256){console.log(`[SFR Auto Detect] Using simple dedupe for ${_.length} candidates`);for(let C=0;C<_.length;C++){const M=_[C];console.log(`[SFR Auto Detect] Simple dedupe candidate ${C+1}/${_.length}`,M.bbox);const S=_.length<=0?1:C/_.length;if(s==null||s(`Detecting candidates: deduplicating (${C}/${_.length})...`,.94+.05*Math.min(1,S)),!k.some(N=>{const E=Math.hypot(M.centerX-N.centerX,M.centerY-N.centerY),X=Math.max(Math.hypot(M.bbox.w,M.bbox.h),Math.hypot(N.bbox.w,N.bbox.h));return Hr(M.bbox,N.bbox)>.28||E<X*.18})&&(k.push(M),k.length>=i))break}return s==null||s("Detecting candidates: deduplicating...",1),k}const v=Math.max(32,Math.round(Math.sqrt(Math.max(1,t*e)/4096))),w=new Map,P=new Set,A=C=>Math.floor(C/v),F=(C,M)=>{if(!Number.isFinite(C.bbox.x)||!Number.isFinite(C.bbox.y)||!Number.isFinite(C.bbox.w)||!Number.isFinite(C.bbox.h)||C.bbox.w<=0||C.bbox.h<=0||C.bbox.w>t*4||C.bbox.h>e*4)return;const S=A(C.bbox.x),I=A(C.bbox.x+C.bbox.w),N=A(C.bbox.y),E=A(C.bbox.y+C.bbox.h);for(let X=N;X<=E;X++)for(let U=S;U<=I;U++){const L=`${U},${X}`,B=w.get(L);B?B.push(M):w.set(L,[M])}};for(let C=0;C<_.length;C++){const M=_[C];if(C===0||C%200===0){const U=_.length<=0?1:C/_.length;s==null||s(`Detecting candidates: deduplicating (${C}/${_.length})...`,.94+.05*Math.min(1,U))}P.clear();const S=A(M.bbox.x),I=A(M.bbox.x+M.bbox.w),N=A(M.bbox.y),E=A(M.bbox.y+M.bbox.h);let X=!1;for(let U=N-1;U<=E+1&&!X;U++)for(let L=S-1;L<=I+1&&!X;L++){const B=w.get(`${L},${U}`);if(B)for(const O of B){if(P.has(O))continue;P.add(O);const z=k[O];if(!z)continue;const J=Math.hypot(M.centerX-z.centerX,M.centerY-z.centerY),V=Math.max(Math.hypot(M.bbox.w,M.bbox.h),Math.hypot(z.bbox.w,z.bbox.h));if(Hr(M.bbox,z.bbox)>.28||J<V*.18){X=!0;break}}}if(!X){const U=k.length;if(k.push(M),F(M,U),k.length>=i)break}}return s==null||s("Detecting candidates: deduplicating...",1),k}function cl(n,t,e){const i=new Uint8Array(n.length);for(let r=0;r<e;r++)for(let s=0;s<t;s++){let o=0,a=0;for(let l=-1;l<=1;l++){const u=r+l;if(!(u<0||u>=e))for(let h=-1;h<=1;h++){const c=s+h;c<0||c>=t||(o+=n[u*t+c],a++)}}i[r*t+s]=Math.round(o/Math.max(1,a))}return i}function ul(n,t,e,i){const r=Math.max(t,e);if(r<=i)return{gray:n,width:t,height:e,scale:1};const s=i/r,o=Math.max(1,Math.round(t*s)),a=Math.max(1,Math.round(e*s)),l=new Uint8Array(o*a);for(let u=0;u<a;u++){const h=Math.min(e-1,Math.floor(u/s));for(let c=0;c<o;c++){const g=Math.min(t-1,Math.floor(c/s));l[u*o+c]=n[h*t+g]}}return{gray:l,width:o,height:a,scale:s}}function Hr(n,t){const e=Math.max(n.x,t.x),i=Math.max(n.y,t.y),r=Math.min(n.x+n.w,t.x+t.w),s=Math.min(n.y+n.h,t.y+t.h),o=Math.max(0,r-e),a=Math.max(0,s-i),l=o*a;if(l<=0)return 0;const u=n.w*n.h+t.w*t.h-l;return u>0?l/u:0}function Qr(n){const t=n.length;if(t===0)return{count:0,mean:0,std:1/0};let e=0;for(const s of n)e+=s;const i=e/t;let r=0;for(const s of n){const o=s-i;r+=o*o}return r/=t,{count:t,mean:i,std:Math.sqrt(Math.max(0,r))}}function hl(n,t,e,i){return{p1:{x:n.x-t*i,y:n.y-e*i},p2:{x:n.x+t*i,y:n.y+e*i}}}function qr(n,t,e,i,r){return[{x:n.p1.x+t*i,y:n.p1.y+e*i},{x:n.p2.x+t*i,y:n.p2.y+e*i},{x:n.p2.x+t*r,y:n.p2.y+e*r},{x:n.p1.x+t*r,y:n.p1.y+e*r}]}function dl(n,t,e){let i=0;for(let r=0;r<4;r++){const s=e[r],o=e[(r+1)%4],a=(o.x-s.x)*(t-s.y)-(o.y-s.y)*(n-s.x);if(Math.abs(a)<=1e-6)continue;const l=a>0?1:-1;if(i===0)i=l;else if(i!==l)return!1}return!0}function Kr(n,t,e,i){const r=Mt(Vt(i,1),t,e);if(!r)return[];const s=[];for(let o=r.y;o<r.y+r.h;o++)for(let a=r.x;a<r.x+r.w;a++)dl(a,o,i)&&s.push(n[o*t+a]);return s}function fl(n,t,e,i,r,s,o,a,l,u,h,c,g,d){const f=h*2,p=c*2,m=Math.hypot(i[1].x-i[0].x,i[1].y-i[0].y),y=Math.hypot(i[2].x-i[1].x,i[2].y-i[1].y),x=Math.hypot(i[2].x-i[3].x,i[2].y-i[3].y),b=Math.hypot(i[3].x-i[0].x,i[3].y-i[0].y),k=Math.max(...[m,y,x,b]),v=Math.max(2,Math.min(f,p)),w=Math.max(4,k*.25),P=Math.max(2,Math.min(12,v*.22)),A=Math.max(1,Math.min(P,Math.max(1,v*.5-1))),F=1,C=Math.max(8,Math.round(Math.min(w,P*3))),M=[[i[0],i[1],i[1],i[0]],[i[1],i[2],i[2],i[1]],[i[2],i[3],i[3],i[2]],[i[3],i[0],i[0],i[3]]],S=[[i[0],i[1],i[1],i[0]],[i[1],i[2],i[2],i[1]],[i[2],i[3],i[3],i[2]],[i[3],i[0],i[0],i[3]]],I=[],N=[],E=[{corners:[i[0],i[1]],seedLine:r,sideLength:m},{corners:[i[1],i[2]],seedLine:s,sideLength:y},{corners:[i[2],i[3]],seedLine:o,sideLength:x},{corners:[i[3],i[0]],seedLine:a,sideLength:b}];for(let G=0;G<E.length;G++){const Y=E[G],et=Math.max(1,Y.sideLength*.5-1),nt=Math.max(1,Math.min(et,w*.5)),it={x:(Y.corners[0].x+Y.corners[1].x)*.5,y:(Y.corners[0].y+Y.corners[1].y)*.5},ct=hl(it,Y.seedLine.dirX,Y.seedLine.dirY,nt),dt=ye(n,t,e,ct.p1,ct.p2,nt,Math.max(4,F+Math.max(P,A)+2)),rt=(dt==null?void 0:dt.line)||ct,Bt=rt.p2.x-rt.p1.x,Xt=rt.p2.y-rt.p1.y,j=Math.hypot(Bt,Xt);if(!Number.isFinite(j)||j<=1e-6)return{ok:!1,score:0,meanSpread:1/0,meanSpreadLimit:1/0,avgStd:1/0,avgStdLimit:1/0,contrast:0,outerMean:0,outerSideMeans:[0,0,0,0],outerSideStds:[1/0,1/0,1/0,1/0],outerSideStdLimit:1/0,outerSideQuads:M,innerSideOk:!1,innerSideStds:[1/0,1/0,1/0,1/0],innerSideStdLimit:1/0,innerSideQuads:S};let ut=-Xt/j,ot=Bt/j;const pt={x:(rt.p1.x+rt.p2.x)*.5,y:(rt.p1.y+rt.p2.y)*.5};(pt.x-l)*ut+(pt.y-u)*ot<0&&(ut=-ut,ot=-ot);const mt=qr(rt,ut,ot,F,F+P),ee=qr(rt,ut,ot,-F,-(F+A));M[G]=mt,S[G]=ee,I.push(Kr(n,t,e,mt)),N.push(Kr(n,t,e,ee))}const X=I.map(Qr);if(X.some(G=>G.count<C||!Number.isFinite(G.std)))return{ok:!1,score:0,meanSpread:1/0,meanSpreadLimit:1/0,avgStd:1/0,avgStdLimit:1/0,contrast:0,outerMean:0,outerSideMeans:[0,0,0,0],outerSideStds:[1/0,1/0,1/0,1/0],outerSideStdLimit:1/0,outerSideQuads:M,innerSideOk:!1,innerSideStds:[1/0,1/0,1/0,1/0],innerSideStdLimit:1/0,innerSideQuads:S};const U=N.map(Qr);if(U.some(G=>G.count<C||!Number.isFinite(G.std)||!Number.isFinite(G.mean)))return{ok:!1,score:0,meanSpread:1/0,meanSpreadLimit:1/0,avgStd:1/0,avgStdLimit:1/0,contrast:0,outerMean:0,outerSideMeans:[0,0,0,0],outerSideStds:[1/0,1/0,1/0,1/0],outerSideStdLimit:1/0,outerSideQuads:M,innerSideOk:!1,innerSideStds:[1/0,1/0,1/0,1/0],innerSideStdLimit:1/0,innerSideQuads:S};const L=X.map(G=>G.mean),B=L.reduce((G,Y)=>G+Y,0)/L.length,O=U.reduce((G,Y)=>G+Y.mean,0)/U.length,z=Math.abs(O-B),J=Math.max(...L)-Math.min(...L),V=X.reduce((G,Y)=>G+Y.std,0)/X.length,q=Math.max(0,d),W=Math.max(6,Math.min(20,z*.45)),tt=L,at=X.map(G=>G.std),D=Math.max(W,Math.min(30,W*g)),Z=U.map(G=>G.std),H=Z.every(G=>G<=D),T=J<=q&&V<=W,Q=1/(1+J/Math.max(1,q)+V/Math.max(1,W));return{ok:T,score:Q,meanSpread:J,meanSpreadLimit:q,avgStd:V,avgStdLimit:W,contrast:z,outerMean:B,outerSideMeans:tt,outerSideStds:at,outerSideStdLimit:W,outerSideQuads:M,innerSideOk:H,innerSideStds:Z,innerSideStdLimit:D,innerSideQuads:S}}function ye(n,t,e,i,r,s,o){const a=r.x-i.x,l=r.y-i.y,u=Math.hypot(a,l);if(!Number.isFinite(u)||u<=1e-6)return null;const h=a/u,c=l/u,g=-c,d=h,f=(i.x+r.x)*.5,p=(i.y+r.y)*.5,m=_e({p1:i,p2:r},o+2);if(!Mt(Vt(m||[i,r],2),t,e))return null;const x=Math.max(8,Math.round(s*2)+1),b=Math.max(8,Math.round(o*2)+1),_=x>1?s*2/(x-1):0,k=b>1?o*2/(b-1):0,v=Array.from({length:x},()=>new Array(b).fill(0));for(let M=0;M<x;M++){const S=-s+_*M;for(let I=0;I<b;I++){const N=-o+k*I,E=f+S*h+N*g,X=p+S*c+N*d;v[M][I]=pl(n,t,e,E,X)}}const w=Ts(v,-o,-s,k,_,!0);if(w.length<8)return null;const P=w.map(M=>{const S=M.x,I=M.y;return{x:f+I*h+S*g,y:p+I*c+S*d,weight:M.weight}}),A=Ce(P);if(!A)return null;let F=A.dirX,C=A.dirY;return F*h+C*c<0&&(F=-F,C=-C),{line:{p1:{x:A.pointX-F*s,y:A.pointY-C*s},p2:{x:A.pointX+F*s,y:A.pointY+C*s}},fitPoints:P.map(M=>({x:M.x,y:M.y}))}}function pl(n,t,e,i,r){if(t<=0||e<=0||n.length!==t*e)return 0;const s=Math.max(0,Math.min(t-1,i)),o=Math.max(0,Math.min(e-1,r)),a=Math.floor(s),l=Math.floor(o),u=Math.min(t-1,a+1),h=Math.min(e-1,l+1),c=s-a,g=o-l,d=n[l*t+a],f=n[l*t+u],p=n[h*t+a],m=n[h*t+u],y=d*(1-c)+f*c,x=p*(1-c)+m*c;return y*(1-g)+x*g}function ml(n,t,e,i,r){if(i<=0||r<=0||i>=t-1||r>=e-1)return{gx:0,gy:0};const s=r*t+i;return{gx:(n[s+1]-n[s-1])*.5,gy:(n[s+t]-n[s-t])*.5}}function gl(n){if(n.length<20)return null;const t=n.map(w=>Math.max(0,w.weight));let e=0;for(const w of t)e=Math.max(e,w);if(!(e>0))return null;for(let w=0;w<t.length;w++)t[w]/=e;const i=w=>{let P=0,A=0,F=0;for(let W=0;W<n.length;W++){const tt=w[W];tt>0&&(P+=tt,A+=n[W].x*tt,F+=n[W].y*tt)}if(!(P>0))return null;A/=P,F/=P;let C=0,M=0,S=0;for(let W=0;W<n.length;W++){const tt=w[W];if(!(tt>0))continue;const at=n[W].x-A,D=n[W].y-F;C+=tt*at*at,M+=tt*at*D,S+=tt*D*D}C/=P,M/=P,S/=P;const I=C+S,N=C*S-M*M,E=-I,X=N,U=Math.max(0,E*E-4*X),L=-.5*(E+(E>=0?1:-1)*Math.sqrt(U)),B=Math.abs(L)>1e-12?L:0,O=Math.abs(L)>1e-12?X/L:I,z=Math.max(B,O);let J=0,V=1;Math.abs(M)>1e-10?(J=z-S,V=M):C>S&&(J=1,V=0);const q=Math.atan2(-J,V);return{centroid:{x:A,y:F},angle:q,totalWeight:P}},r=i(t);if(!r)return null;const s=Math.cos(r.angle),o=Math.sin(r.angle),a=new Array(2*16*8).fill(0);for(let w=0;w<n.length;w++){const P=n[w].x-r.centroid.x,A=n[w].y-r.centroid.y,F=P*s+A*o,C=Math.round(F*8+16*8);if(C>=3&&C<a.length-3)for(let M=-3;M<=3;M++)a[C+M]+=t[w]}let l=16*8;for(let w=-5*8+16*8;w<=5*8+16*8;w++)a[w]>a[l]&&(l=w);let u=l-1;for(;u>1&&a[u]>.05*a[l];)u--;let h=l+1;for(;h<a.length-1&&a[h]>.05*a[l];)h++;let c=Math.max(1,u-8);for(;c>1&&a[c]<=a[u];)c--;let g=Math.min(a.length-1,h+8);for(;g<a.length-1&&a[g]<=a[h];)g++;const d=a.slice();for(let w=1;w<d.length;w++)d[w]+=d[w-1];const f=d[d.length-1];if(!(f>0))return null;let p=0;for(let w=1;w<d.length;w++)Math.abs(d[w]-.1*f)<Math.abs(d[p]-.1*f)&&(p=w);let m=d.length-1;for(let w=d.length-2;w>0;w--)Math.abs(d[w]-.9*f)<Math.abs(d[m]-.9*f)&&(m=w);let y=p/8-16,x=m/8-16;const b=x-y;y-=b*.7,x+=b*.7,y=Math.max((c+u)/16-16,y),x=Math.min((g+h)/16-16,x);const _=t.slice();for(let w=0;w<n.length;w++){const P=n[w].x-r.centroid.x,A=n[w].y-r.centroid.y,F=P*s+A*o;_[w]=F>=y&&F<=x?t[w]**4*(1/(10+Math.abs(F))):0}const k=i(_);if(!k)return null;const v=[];for(let w=0;w<n.length;w++)_[w]>0&&v.push({x:n[w].x,y:n[w].y,weight:_[w]});return v.length<8?null:{centroid:k.centroid,angle:k.angle,keptSamples:v}}function yl(n,t,e,i,r,s=Lt){var P;const o=r.x-i.x,a=r.y-i.y,l=Math.hypot(o,a);if(!Number.isFinite(l)||l<=12)return null;const u=o/l,h=a/l,c=5,g=4*s+.5,d=(A,F,C,M,S)=>{const I={x:A.x-F*l*.5,y:A.y-C*l*.5},N={p1:I,p2:{x:I.x+F*l,y:I.y+C*l}},E=_e(N,g+2),X=Mt(Vt(E??[N.p1,N.p2],3),t,e),U=[],L=new Map;if(!X)return{reduced:null,scanlines:L};for(let B=X.y;B<X.y+X.h;B++)for(let O=X.x;O<X.x+X.w;O++){const z=O,J=B,V=z-I.x,q=J-I.y,W=V*F+q*C;if(!(W>c&&W<l-c))continue;const tt=z-A.x,at=J-A.y,D=tt*M+at*S;if(Math.abs(D)<12){const{gx:Z,gy:H}=ml(n,t,e,O,B),T=Z*Z+H*H;T>0&&U.push({x:z,y:J,weight:T})}if(Math.abs(D)<g){const Z=L.get(B);Z?(O<Z.start&&(Z.start=O),O>Z.end&&(Z.end=O)):L.set(B,{start:O,end:O})}}return{reduced:gl(U),scanlines:L}};let f={x:(i.x+r.x)*.5,y:(i.y+r.y)*.5},p=u,m=h,y=-m,x=p,b=d(f,p,m,y,x);if(!b.reduced)return null;f=b.reduced.centroid,y=Math.cos(b.reduced.angle),x=Math.sin(b.reduced.angle),p=-x,m=y,p*u+m*h<0&&(p=-p,m=-m,y=-y,x=-x);let _=d(f,p,m,y,x);if(!_.reduced)return null;const k=Math.hypot(_.reduced.centroid.x-f.x,_.reduced.centroid.y-f.y);f=_.reduced.centroid,y=Math.cos(_.reduced.angle),x=Math.sin(_.reduced.angle),p=-x,m=y,p*u+m*h<0&&(p=-p,m=-m,y=-y,x=-x);let v=_;if(k>1){const A=d(f,p,m,y,x);A.reduced&&(v=A,f=A.reduced.centroid,y=Math.cos(A.reduced.angle),x=Math.sin(A.reduced.angle),p=-x,m=y,p*u+m*h<0&&(p=-p,m=-m))}const w=(((P=v.reduced)==null?void 0:P.keptSamples)??[]).map(A=>({x:A.x,y:A.y}));return w.length<8?null:{line:{p1:{x:f.x-p*l*.5,y:f.y-m*l*.5},p2:{x:f.x+p*l*.5,y:f.y+m*l*.5}},fitPoints:w,correctedScanlines:v.scanlines}}function xl(n,t,e){var b;const i=n.length,r=((b=n[0])==null?void 0:b.length)??0;if(r===0||i===0)return 0;const s=Math.max(0,Math.min(r-1,t)),o=Math.max(0,Math.min(i-1,e)),a=Math.floor(s),l=Math.floor(o),u=Math.min(r-1,a+1),h=Math.min(i-1,l+1),c=s-a,g=o-l,d=n[l][a],f=n[l][u],p=n[h][a],m=n[h][u],y=d*(1-c)+f*c,x=p*(1-c)+m*c;return y*(1-g)+x*g}function bl(n,t,e,i,r,s,o){var U;const a=i.p2.x-i.p1.x,l=i.p2.y-i.p1.y,u=Math.hypot(a,l);if(!Number.isFinite(u)||u<=1e-6)return null;const h=a/u,c=l/u,g=-c,d=h,f=(i.p1.x+i.p2.x)*.5,p=(i.p1.y+i.p2.y)*.5,m=_e(i,s+2),y=Mt(Vt(m||[i.p1,i.p2],3),t,e);if(!y)return null;const x=Ao(n,t,e,y,0,0,dn(o==null?void 0:o.bayerPattern,"RAW edge refinement"),o==null?void 0:o.greenPhase,o==null?void 0:o.blackLevel),b=x.length;if((((U=x[0])==null?void 0:U.length)??0)<6||b<6)return null;const k=Math.max(8,Math.round(r*2)+1),v=Math.max(8,Math.round(s*2)+1),w=k>1?r*2/(k-1):0,P=v>1?s*2/(v-1):0,A=Array.from({length:k},()=>new Array(v).fill(0));for(let L=0;L<k;L++){const B=-r+w*L;for(let O=0;O<v;O++){const z=-s+P*O,J=f+B*h+z*g,V=p+B*c+z*d;A[L][O]=xl(x,J-y.x,V-y.y)}}const{gx:F,gy:C}=To(A),M=F>=C,S=Ts(A,-s,-r,P,w,M);if(S.length<8)return null;const I=S.map(L=>{const B=L.x,O=L.y;return{x:f+O*h+B*g,y:p+O*c+B*d,weight:L.weight}}),N=Ce(I);if(!N)return null;let E=N.dirX,X=N.dirY;return E*h+X*c<0&&(E=-E,X=-X),{line:{p1:{x:N.pointX-E*r,y:N.pointY-X*r},p2:{x:N.pointX+E*r,y:N.pointY+X*r}},fitPoints:I.map(L=>({x:L.x,y:L.y}))}}function ke(n){const t=Math.max(0,Math.min(1,n));return t<=.04045?t/12.92:Math.pow((t+.055)/1.055,2.4)}function or(n,t,e){if(t<=0||e<=0||n.length!==t*e)return new Uint8Array(Math.max(0,t*e));let i=1/0,r=-1/0;for(let f=0;f<n.length;f++){const p=n[f];Number.isFinite(p)&&(p<i&&(i=p),p>r&&(r=p))}if(!Number.isFinite(i)||!Number.isFinite(r)||r<=i+1e-9)return new Uint8Array(t*e);const s=1024,o=new Uint32Array(s),a=r-i;for(let f=0;f<n.length;f++){const p=Math.max(0,Math.min(1,(n[f]-i)/a)),m=Math.min(s-1,Math.max(0,Math.floor(p*(s-1))));o[m]++}const l=n.length,u=f=>{const p=l*f;let m=0;for(let y=0;y<s;y++)if(m+=o[y],m>=p)return i+y/Math.max(1,s-1)*a;return r},h=u(.01),c=u(.99),g=Math.max(1e-9,c-h),d=new Uint8Array(t*e);for(let f=0;f<n.length;f++){const p=Math.max(0,Math.min(1,(n[f]-h)/g));d[f]=Math.round(p*255)}return d}function $n(n,t,e=0){const i=new Float32Array(n.width*n.height),r=n.data;for(let s=0,o=0;s<r.length;s+=4,o++)i[o]=Ls(r,s,t,e);return or(i,n.width,n.height)}function _l(n){return Number.isFinite(n)?Math.max(0,Math.min(65535,Number(n))):0}function Ls(n,t,e,i=0){let r=n[t]/255,s=n[t+1]/255,o=n[t+2]/255;e&&(r=ke(r),s=ke(s),o=ke(o));const a=.2126*r+.7152*s+.0722*o;return Math.max(0,a-_l(i)/65535)}function Ml(n,t){const e=n.width,i=n.height,r=n.data;if(r.length<e*i*3)return new Uint8Array(e*i);const s=new Float32Array(e*i);for(let o=0;o<e*i;o++){const a=o*3;t!==void 0?s[o]=r[a+t]:s[o]=.2126*r[a]+.7152*r[a+1]+.0722*r[a+2]}return or(s,e,i)}function wl(n){const t=new Float32Array(n.width*n.height);for(let e=0;e<n.data.length;e++)t[e]=n.data[e];return or(t,n.width,n.height)}function $r(n,t,e){const i=Mt(t,n.width,n.height);if(!i)return null;const r=new Uint16Array(i.w*i.h*3),s=n.data;let o=0;for(let a=i.y;a<i.y+i.h;a++)for(let l=i.x;l<i.x+i.w;l++){const u=(a*n.width+l)*4;let h=s[u]/255,c=s[u+1]/255,g=s[u+2]/255;e&&(h=ke(h),c=ke(c),g=ke(g)),r[o++]=Math.max(0,Math.min(65535,Math.round(h*65535))),r[o++]=Math.max(0,Math.min(65535,Math.round(c*65535))),r[o++]=Math.max(0,Math.min(65535,Math.round(g*65535)))}return{data:r,width:i.w,height:i.h}}function Jr(n,t,e,i=0){const r=Mt(t,n.width,n.height);if(!r)return null;const s=new Uint16Array(r.w*r.h*3),o=n.data;let a=0;for(let l=r.y;l<r.y+r.h;l++)for(let u=r.x;u<r.x+r.w;u++){const h=(l*n.width+u)*4,c=Math.max(0,Math.min(65535,Math.round(Ls(o,h,e,i)*65535)));s[a++]=c,s[a++]=c,s[a++]=c}return{data:s,width:r.w,height:r.h}}function Sl(n,t){const e=Mt(t,n.width,n.height);if(!e)return null;const i=new Uint16Array(e.w*e.h),r=n.data;let s=0;for(let o=e.y;o<e.y+e.h;o++)for(let a=e.x;a<e.x+e.w;a++){const l=(o*n.width+a)*4;i[s++]=Math.max(0,Math.min(65535,Math.round((.2126*r[l]+.7152*r[l+1]+.0722*r[l+2])*257)))}return{data:i,width:e.w,height:e.h}}function vl(n,t){const e=Mt(t,n.width,n.height);if(!e)return null;const i=new Uint16Array(e.w*e.h);let r=0;for(let s=e.y;s<e.y+e.h;s++){const o=s*n.width;for(let a=e.x;a<e.x+e.w;a++)i[r++]=n.data[o+a]}return{data:i,width:e.w,height:e.h}}function qt(n,t,e){return{x:n.x*t,y:n.y*e}}function Pl(n,t,e){return{p1:qt(n.p1,t,e),p2:qt(n.p2,t,e)}}function Hi(n,t){const e=t(n);return{x:Number.isFinite(e.x)?e.x:n.x,y:Number.isFinite(e.y)?e.y:n.y}}function Cl(n,t){return n.map(e=>Hi(e,t))}function ge(n,t,e,i=0,r=0){if(!n||n.length<8)return;const s=n.map(o=>({x:o.x*t-i,y:o.y*e-r})).filter(o=>Number.isFinite(o.x)&&Number.isFinite(o.y));return s.length>=8?s:void 0}function kl(n,t){return{p1:Hi(n.p1,t),p2:Hi(n.p2,t)}}function fn(n,t,e){return{p1:{x:n.p1.x-t,y:n.p1.y-e},p2:{x:n.p2.x-t,y:n.p2.y-e}}}function Fl(n,t,e,i){const r=Math.max(0,Math.min(n.width-1,t)),o=(Math.max(0,Math.min(n.height-1,e))*n.width+r)*4;let a=n.data[o]/255,l=n.data[o+1]/255,u=n.data[o+2]/255;return i&&(a=ke(a),l=ke(l),u=ke(u)),(.2126*a+.7152*l+.0722*u)*65535}function Es(n){return n.kind==="u16-mono"}function sn(n){return n.width}function an(n){return n.height}function ai(n,t,e,i){if(Es(n)){const r=Math.max(0,Math.min(n.width-1,t)),s=Math.max(0,Math.min(n.height-1,e));return n.data[s*n.width+r]}return Fl(n,t,e,i)}function Al(n,t,e,i){if(Es(n)&&n.coordinateSpace==="distorted-padded"){const r=Math.round(n.paddingOffsetX??0),s=Math.round(n.paddingOffsetY??0);return ai(n,t+r,e+s,i)}return ai(n,t,e,i)}function Tl(n,t,e,i,r=3){const o=[...n,{x:(n[0].x+n[1].x+n[2].x+n[3].x)*.25,y:(n[0].y+n[1].y+n[2].y+n[3].y)*.25},{x:(n[0].x+n[1].x)*.5,y:(n[0].y+n[1].y)*.5},{x:(n[1].x+n[2].x)*.5,y:(n[1].y+n[2].y)*.5},{x:(n[2].x+n[3].x)*.5,y:(n[2].y+n[3].y)*.5},{x:(n[3].x+n[0].x)*.5,y:(n[3].y+n[0].y)*.5}].map(a=>Fe(a,t)).filter(a=>Number.isFinite(a.x)&&Number.isFinite(a.y));return o.length===0?null:Mt(Vt(o,r),e,i)}function lr(n,t,e,i){const r=new Map;for(let s=n.y;s<n.y+n.h;s++)for(let o=n.x;o<n.x+n.w;o++){const a=Ut({x:o,y:s},t);if(!Number.isFinite(a.x)||!Number.isFinite(a.y))continue;const l=Math.round(a.x),u=Math.round(a.y);if(l<0||u<0||l>=e||u>=i)continue;const h=r.get(u);h?(l<h.start&&(h.start=l),l>h.end&&(h.end=l)):r.set(u,{start:l,end:l})}return r}function Il(n,t,e,i,r,s,o){const a=new Map,l=t.p2.x-t.p1.x,u=t.p2.y-t.p1.y,h=Math.hypot(l,u);if(!Number.isFinite(h)||h<=1e-6)return a;const c=l/h,g=u/h,d=-g,f=c,p={x:(t.p1.x+t.p2.x)*.5,y:(t.p1.y+t.p2.y)*.5},m=Math.max(1,e+1),y=Math.max(1,i+1.5);for(let x=n.y;x<n.y+n.h;x++)for(let b=n.x;b<n.x+n.w;b++){const _=b+.5,k=x+.5,v=_-p.x,w=k-p.y,P=v*c+w*g;if(!Number.isFinite(P)||Math.abs(P)>m)continue;const A=v*d+w*f;if(!Number.isFinite(A)||Math.abs(A)>y)continue;const F=Ut({x:_,y:k},r);if(!Number.isFinite(F.x)||!Number.isFinite(F.y))continue;const C=Math.round(F.x),M=Math.round(F.y);if(C<0||M<0||C>=s||M>=o)continue;const S=a.get(M);S?(C<S.start&&(S.start=C),C>S.end&&(S.end=C)):a.set(M,{start:C,end:C})}return a}function Us(n,t,e,i,r,s){const o=new Map,a=t.p2.x-t.p1.x,l=t.p2.y-t.p1.y,u=Math.hypot(a,l);if(!Number.isFinite(u)||u<=1e-6)return o;const h=a/u,c=l/u,g=-c,d=h,f={x:(t.p1.x+t.p2.x)*.5,y:(t.p1.y+t.p2.y)*.5},p=Math.max(1,e+1),m=Math.max(1,i+1.5),y=Mt(n,r,s);if(!y)return o;for(let x=y.y;x<y.y+y.h;x++)for(let b=y.x;b<y.x+y.w;b++){const _=b-f.x,k=x-f.y,v=_*h+k*c;if(!Number.isFinite(v)||Math.abs(v)>p)continue;const w=_*g+k*d;if(!Number.isFinite(w)||Math.abs(w)>m)continue;const P=o.get(x);P?(b<P.start&&(P.start=b),b>P.end&&(P.end=b)):o.set(x,{start:b,end:b})}return o}function Bs(n,t,e,i){const r=new Map;for(const[s,o]of n)for(let a=o.start;a<=o.end;a++){const l=Fe({x:a,y:s},t);if(!Number.isFinite(l.x)||!Number.isFinite(l.y))continue;const u=Math.round(l.x),h=Math.round(l.y);if(u<0||h<0||u>=e||h>=i)continue;const c=r.get(h);c?(u<c.start&&(c.start=u),u>c.end&&(c.end=u)):r.set(h,{start:u,end:u})}return r}function cr(n){return Math.abs(n.k1)<1e-4&&Math.abs(n.k2)<1e-4}function Rl(n){return[{x:n.x,y:n.y},{x:n.x+n.w,y:n.y},{x:n.x+n.w,y:n.y+n.h},{x:n.x,y:n.y+n.h}]}function Kt(n,t,e,i,r){return Fe({x:i.x+n*t,y:i.y+n*e},r)}function ur(n,t,e,i,r){const o=Kt(n,t,e,i,r),a=Kt(n+1e-4,t,e,i,r);return{x:(a.x-o.x)/1e-4,y:(a.y-o.y)/1e-4}}function Ds(n,t,e,i,r,s){let o=.01;const a=h=>{const c=Kt(h,t,e,i,s);return Math.hypot(c.x-r.x,c.y-r.y)},l=a(n),u=a(n+o);if(!Number.isFinite(l)||!Number.isFinite(u))return null;if(l>u){let h=n,c=n+o;for(let g=0;g<24;g++){o*=2;const d=h+o,f=a(d),p=a(c);if(!Number.isFinite(f)||!Number.isFinite(p))break;if(f>=p)return{a:h,b:d};h=c,c=d}}else{let h=n,c=n+o;for(let g=0;g<24;g++){o*=2;const d=c-o,f=a(d),p=a(h);if(!Number.isFinite(f)||!Number.isFinite(p))break;if(f>=p)return{a:d,b:c};c=h,h=d}}return{a:n-Math.max(.5,o),b:n+Math.max(.5,o)}}function Nl(n,t,e=33){const i=n.p2.x-n.p1.x,r=n.p2.y-n.p1.y,s=Math.hypot(i,r);if(!Number.isFinite(s)||s<=1e-6)return[Fe(n.p1,t),Fe(n.p2,t)];const o=i/s,a=r/s,l={x:(n.p1.x+n.p2.x)*.5,y:(n.p1.y+n.p2.y)*.5},u=s*.5,h=Math.max(9,e),c=[];for(let g=0;g<h;g++){const d=h===1?.5:g/(h-1),f=-u+d*(u*2);c.push(Kt(f,o,a,l,t))}return c}function Ll(n,t,e,i,r,s,o=1){const a=n.p2.x-n.p1.x,l=n.p2.y-n.p1.y,u=Math.hypot(a,l);if(!Number.isFinite(u)||u<=1e-6)return null;const h=a/u,c=l/u,g={x:(n.p1.x+n.p2.x)*.5,y:(n.p1.y+n.p2.y)*.5},d=Math.max(24,Math.round(e*2)+1),f=[];for(let p=0;p<d;p++){const m=d===1?.5:p/(d-1),y=-e+m*(e*2),x=Kt(y,h,c,g,t),b=ur(y,h,c,g,t),_=Math.hypot(b.x,b.y);if(!Number.isFinite(_)||_<=1e-9)continue;const k=-b.y/_,v=b.x/_;f.push({x:x.x+k*(i+o),y:x.y+v*(i+o)},{x:x.x-k*(i+o),y:x.y-v*(i+o)})}if(f.length<2){const p={p1:Fe(n.p1,t),p2:Fe(n.p2,t)},m=_e(p,i+o);return m?Mt(Vt(m,o),r,s):null}return Mt(Vt(f,o),r,s)}function rn(n,t,e,i,r){return Ut({x:i.x+n*t,y:i.y+n*e},r)}function El(n,t,e,i,r,s){let o=.01;const a=h=>{const c=rn(h,t,e,i,s);return Math.hypot(c.x-r.x,c.y-r.y)},l=a(n),u=a(n+o);if(!Number.isFinite(l)||!Number.isFinite(u))return null;if(l>u){let h=n,c=n+o;for(let g=0;g<24;g++){o*=2;const d=h+o,f=a(d),p=a(c);if(!Number.isFinite(f)||!Number.isFinite(p))break;if(f>=p)return{a:h,b:d};h=c,c=d}}else{let h=n,c=n+o;for(let g=0;g<24;g++){o*=2;const d=c-o,f=a(d),p=a(h);if(!Number.isFinite(f)||!Number.isFinite(p))break;if(f>=p)return{a:d,b:c};c=h,h=d}}return{a:n-Math.max(.5,o),b:n+Math.max(.5,o)}}function hr(n,t,e){const i=(t.x-n.x)*(t.y-e.y)-(t.x-e.x)*(t.y-n.y);if(Math.abs(i)<=1e-12)return .5*(n.x+e.x);const r=(t.x-n.x)*(t.x-n.x)*(t.y-e.y)-(t.x-e.x)*(t.x-e.x)*(t.y-n.y),s=t.x-.5*r/i;return Number.isFinite(s)?s:.5*(n.x+e.x)}function Ul(n,t,e,i,r){const o=rn(n,t,e,i,r),a=rn(n+1e-4,t,e,i,r);return{x:(a.x-o.x)/1e-4,y:(a.y-o.y)/1e-4}}function Bl(n,t,e,i,r,s,o=!1,a){const l=[],u=[],h=a?Lt*2:Lt,c=Math.max(1,Math.min(s,h)),g=i.p2.x-i.p1.x,d=i.p2.y-i.p1.y,f=Math.hypot(g,d);if(!Number.isFinite(f)||f<=1e-6)return null;const p=g/f,m=d/f,y={x:(i.p1.x+i.p2.x)*.5,y:(i.p1.y+i.p2.y)*.5},x={p1:Ut(i.p1,e),p2:Ut(i.p2,e)},b=x.p2.x-x.p1.x,_=x.p2.y-x.p1.y,k=Math.hypot(b,_);if(!Number.isFinite(k)||k<=1e-6)return null;const v=b/k,w=_/k,P=-w,A=v,F={x:(x.p1.x+x.p2.x)*.5,y:(x.p1.y+x.p2.y)*.5},C=Mt(a||Vt(_e(i,s+2)??[i.p1,i.p2],2),n.width,n.height);if(!C)return null;const M=Il(C,i,r,c,e,sn(t),an(t));if(M.size===0)return null;const S=!cr(e);for(const[L,B]of M)if(!(L<0||L>=an(t)))for(let O=B.start;O<=B.end;O++){if(O<0||O>=sn(t))continue;const z={x:O,y:L};let J,V;if(S){const q=Fe(z,e);if(!Number.isFinite(q.x)||!Number.isFinite(q.y)||Math.round(q.x)<0||Math.round(q.x)>=n.width||Math.round(q.y)<0||Math.round(q.y)>=n.height)continue;const W=q.x-y.x,tt=q.y-y.y,at=W*p+tt*m;if(!Number.isFinite(at))continue;J=at,V=W*-m+tt*p;const D=El(at,p,m,y,z,e);if(!D)continue;const Z=.5*(D.a+D.b),H=rn(D.a,p,m,y,e),T=rn(Z,p,m,y,e),Q=rn(D.b,p,m,y,e),G=hr({x:D.a,y:Math.hypot(H.x-z.x,H.y-z.y)},{x:Z,y:Math.hypot(T.x-z.x,T.y-z.y)},{x:D.b,y:Math.hypot(Q.x-z.x,Q.y-z.y)});if(!Number.isFinite(G))continue;J=G;const Y=Ul(G,p,m,y,e),et=Math.hypot(Y.x,Y.y);if(!Number.isFinite(et)||et<=1e-9)continue;const nt=Y.x/et,ct=-(Y.y/et),dt=nt,rt=rn(G,p,m,y,e);V=(z.x-rt.x)*ct+(z.y-rt.y)*dt}else{const q=z.x-F.x,W=z.y-F.y;J=q*v+W*w,V=q*P+W*A}!Number.isFinite(J)||Math.abs(J)>r||!Number.isFinite(V)||Math.abs(V)>c||(l.push(V),u.push(ai(t,O,L,o)))}if(l.length<8)return null;const I=Ut(i.p1,e),N=Ut(i.p2,e),E=N.x-I.x,X=N.y-I.y,U=Math.abs(E)>=Math.abs(X)?1:2;return wn(l,u,U,h)}function Dl(n,t,e,i,r,s,o=!1,a,l,u,h=!1){const c=[],g=[],d=a?Lt*2:Lt,f=Math.max(1,Math.min(s,d)),p=i.p2.x-i.p1.x,m=i.p2.y-i.p1.y,y=Math.hypot(p,m);if(!Number.isFinite(y)||y<=1e-6)return null;const x=p/y,b=m/y,_=-b,k=x,v={x:(i.p1.x+i.p2.x)*.5,y:(i.p1.y+i.p2.y)*.5},w=Mt(a||Vt(_e(i,d*4+2)??[i.p1,i.p2],2),n.width,n.height);if(!w)return null;const P=l??(u?lr(Mt(u,sn(t),an(t))??u,e,n.width,n.height):Us(w,i,Math.max(1,r),f*4+.5,n.width,n.height));if(P.size===0)return null;const A=Bs(P,e,sn(t),an(t));if(A.size===0)return null;const F=!cr(e);for(const[M,S]of A)for(let I=S.start;I<=S.end;I++){const N={x:I,y:M},E=Ut(N,e);if(!Number.isFinite(E.x)||!Number.isFinite(E.y)||Math.round(E.x)<0||Math.round(E.x)>=n.width||Math.round(E.y)<0||Math.round(E.y)>=n.height)continue;const X=E.x-v.x,U=E.y-v.y,L=X*x+U*b;let B=X*_+U*k;if(F){const O=Ds(L,x,b,v,N,e);if(!O)continue;const z=.5*(O.a+O.b),J=Kt(O.a,x,b,v,e),V=Kt(z,x,b,v,e),q=Kt(O.b,x,b,v,e),W=hr({x:O.a,y:Math.hypot(J.x-N.x,J.y-N.y)},{x:z,y:Math.hypot(V.x-N.x,V.y-N.y)},{x:O.b,y:Math.hypot(q.x-N.x,q.y-N.y)});if(!Number.isFinite(W))continue;const tt=ur(W,x,b,v,e),at=Math.hypot(tt.x,tt.y);if(!Number.isFinite(at)||at<=1e-9)continue;const D=tt.x/at,H=-(tt.y/at),T=D,Q=Kt(W,x,b,v,e);B=(N.x-Q.x)*H+(N.y-Q.y)*T}!Number.isFinite(L)||Math.abs(L)>Math.max(1,r)||!Number.isFinite(B)||Math.abs(B)>f||(c.push(B),g.push(Al(t,I,M,o)))}if(c.length<8)return null;const C=Math.abs(p)>=Math.abs(m)?1:2;return h?Rn(c,g,C,d):wn(c,g,C,d)}function Ol(n,t,e,i,r,s){const o=n.width,a=n.height,l=dn(n.bayerPattern,"corrected RAW edge measurement"),u=s!=null&&s.correctedRect?Lt*2:Lt,h=Math.max(1,Math.min(r,u)),c=(s==null?void 0:s.restrictToStrip)??!0,g=e.p2.x-e.p1.x,d=e.p2.y-e.p1.y,f=Math.hypot(g,d);if(!Number.isFinite(f)||f<=1e-6)return null;const p=g/f,m=d/f,y=-m,x=p,b={x:(e.p1.x+e.p2.x)*.5,y:(e.p1.y+e.p2.y)*.5},_={p1:{x:b.x-p*Math.max(1,i),y:b.y-m*Math.max(1,i)},p2:{x:b.x+p*Math.max(1,i),y:b.y+m*Math.max(1,i)}},k=_e(_,h+2),v=(s!=null&&s.fixedRawRect?Mt(s.fixedRawRect,o,a):null)??(s!=null&&s.correctedRect?sr(s.correctedRect,t,o,a):null)??(k?Tl(k,t,o,a,2):null);if(!v)return null;const w=[],P=[];for(let C=v.y;C<v.y+v.h;C++){const M=C*o;for(let S=v.x;S<v.x+v.w;S++){if(!vt(S,C,l,s==null?void 0:s.greenPhase))continue;const I=Ut({x:S,y:C},t);if(!Number.isFinite(I.x)||!Number.isFinite(I.y))continue;const N=I.x-b.x,E=I.y-b.y,X=N*p+E*m;if(!Number.isFinite(X)||c&&Math.abs(X)>Math.max(1,i))continue;const U=N*y+E*x;if(!Number.isFinite(U)||c&&Math.abs(U)>h)continue;w.push(U);let L;L=Math.max(0,n.data[M+S]-Mn(s==null?void 0:s.blackLevel,S,C)),P.push(L)}}if(w.length<8)return null;const A=Math.abs(g)>=Math.abs(d)?1:2,F=Math.max(2,(s==null?void 0:s.shortSidePxOverride)??(c?h*2:Math.min(v.w,v.h)));return pi(w,P,F,s==null?void 0:s.manualBinSize,A,s==null?void 0:s.preferAutoPerEdgeBin,!1,!!(s!=null&&s.forceLegacyModel))}function Os(n,t,e,i,r,s=!1,o,a,l,u=!1){if(a&&l){const P=l.p2.x-l.p1.x,A=l.p2.y-l.p1.y,F=Math.hypot(P,A);if(Number.isFinite(F)&&F>1e-6)return Bl(a,n,t,l,Math.max(1,F*.5),r,s,o)}const h=[],c=[],g=Lt,d=e.p2.x-e.p1.x,f=e.p2.y-e.p1.y,p=Math.hypot(d,f);if(!Number.isFinite(p)||p<=1e-6)return null;const m=d/p,y=f/p,x={x:(e.p1.x+e.p2.x)*.5,y:(e.p1.y+e.p2.y)*.5},b=(o?Mt(o,sn(n),an(n)):null)??Ll(e,t,i+1,r+1,sn(n),an(n),1);if(!b)return null;const _=lr(b,t,sn(n),an(n));if(_.size===0)return null;const k=-y,v=m;for(const[P,A]of _)for(let F=A.start;F<=A.end;F++){const C={x:F,y:P};let M=(C.x-x.x)*m+(C.y-x.y)*y,S=(C.x-x.x)*k+(C.y-x.y)*v;!Number.isFinite(M)||Math.abs(M)>i+1||!Number.isFinite(S)||Math.abs(S)>=g||(h.push(S),c.push(ai(n,F,P,s)))}if(h.length<8)return null;const w=Math.abs(d)>=Math.abs(f)?1:2;return u?Rn(h,c,w,g):wn(h,c,w,g)}function Vl(n,t,e){var g;const i=e.sourceMode??(t.isThreePlane?"three-plane":"rggb-raw"),r=e.useQuadraticProjection!==!1,s=!!e.forceRenderedMeasurement,o=n.width,a=n.height,l=e.threePlaneChannel,u=ar(e.detectionTuning),h=e.monochromeBlackLevel??0;if(i==="rggb-raw"&&!s){if(!t||t.isThreePlane)return null;const d=Qo(t,e.greenPhase),f=o/Math.max(1,t.width),p=a/Math.max(1,t.height);return{sourceMode:i,detectionGray:d,detectionWidth:t.width,detectionHeight:t.height,detectToMeasureX:1,detectToMeasureY:1,detectToDisplayX:f,detectToDisplayY:p,measureToDisplayX:f,measureToDisplayY:p,detectPointToDisplay:m=>qt(m,f,p),measurePointToDisplay:m=>qt(m,f,p),displayPointToDetect:m=>qt(m,1/Math.max(1e-9,f),1/Math.max(1e-9,p)),measureUsesDisplayLine:!1,measureWidth:t.width,measureHeight:t.height,refineLine:(m,y,x)=>bl(t.data,t.width,t.height,{p1:m,p2:y},x*.5,Math.max(4,x*.2),{greenPhase:e.greenPhase,bayerPattern:t.bayerPattern})||ye(d,t.width,t.height,m,y,x*.5,Math.max(4,x*.2)),measureEdge:(m,y,x,b,_)=>Pe(t.data,t.width,t.height,m,y,x,b,{greenOnly:!0,greenPhase:e.greenPhase,bayerPattern:t.bayerPattern,blackLevel:e.blackLevel??void 0,preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?ge(_==null?void 0:_.fitPoints,1,1):void 0})}}if(s){const d=!!e.distortionCurveApplied&&!!e.distortionModel,f=i==="rggb-raw"&&!!e.distortionCorrected&&!!e.distortionModel&&!t.isThreePlane,p=!!e.distortionCorrected&&!!e.distortionModel&&!!e.distortionOriginalSamplingPlane,m=!!e.distortionCorrected&&!!e.distortionSamplingPlane,y=n,x=$n(y,!!e.sfrHasGamma,i==="unmix-bw"?h:0);return{sourceMode:i,detectionGray:x,detectionWidth:y.width,detectionHeight:y.height,detectToMeasureX:1,detectToMeasureY:1,detectToDisplayX:1,detectToDisplayY:1,measureToDisplayX:1,measureToDisplayY:1,detectPointToDisplay:b=>b,measurePointToDisplay:b=>b,displayPointToDetect:b=>b,measureUsesDisplayLine:!1,measureWidth:y.width,measureHeight:y.height,refineLine:(b,_,k)=>(p?yl(x,n.width,n.height,b,_,Lt):null)||ye(x,n.width,n.height,b,_,k*.5,Math.max(4,k*.2)),measureEdge:(b,_,k,v,w)=>{const P=e.distortionModel?sr(b,e.distortionModel,t.width,t.height):null;if(f){const M={p1:Ut(_.p1,e.distortionModel),p2:Ut(_.p2,e.distortionModel)},S=Math.hypot(M.p2.x-M.p1.x,M.p2.y-M.p1.y),I=Math.max(2,S*.5*u.sampleHalfWidthRatio);return ho(t,e.distortionModel,M,Math.max(1,S*.5),I,{greenPhase:e.greenPhase,blackLevel:e.blackLevel??void 0,correctedRect:b})}if(f)return Ol(t,e.distortionModel,_,k,v,{greenPhase:e.greenPhase,blackLevel:e.blackLevel??void 0,correctedRect:b,fixedRawRect:P,preferAutoPerEdgeBin:!0});if(p)return Dl(n,e.distortionOriginalSamplingPlane,e.distortionModel,_,k,v,!!e.sfrHasGamma,b,(w==null?void 0:w.correctedScanlines)??null,P);if(d){const M={p1:Ut(_.p1,e.distortionModel),p2:Ut(_.p2,e.distortionModel)},S=Math.hypot(M.p2.x-M.p1.x,M.p2.y-M.p1.y);return Os(e.distortionOriginalSamplingPlane??e.distortionSamplingPlane??e.distortionSamplingImage??n,e.distortionModel,M,Math.max(1,S*.5),v,!!e.sfrHasGamma,b,e.distortionBaseImage??n,_)}if(m){const M=vl(e.distortionSamplingPlane,b);if(!M)return null;const S=fn(_,b.x,b.y);return Pe(M.data,M.width,M.height,{x:0,y:0,w:M.width,h:M.height},S,k,v,{preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?ge(w==null?void 0:w.fitPoints,1,1,b.x,b.y):void 0})}const F=i==="unmix-bw"?Jr(n,b,!!e.sfrHasGamma,h):$r(n,b,!!e.sfrHasGamma);if(!F)return null;const C=fn(_,b.x,b.y);return Pe(F.data,F.width,F.height,{x:0,y:0,w:F.width,h:F.height},C,k,v,{isThreePlane:!0,threePlaneChannel:void 0,preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?ge(w==null?void 0:w.fitPoints,1,1,b.x,b.y):void 0})}}}if(i==="three-plane"){if(t.isThreePlane&&!e.sfrHasGamma){const f=Ml(t,l),p=o/Math.max(1,t.width),m=a/Math.max(1,t.height);return{sourceMode:i,detectionGray:f,detectionWidth:t.width,detectionHeight:t.height,detectToMeasureX:1,detectToMeasureY:1,detectToDisplayX:p,detectToDisplayY:m,measureToDisplayX:p,measureToDisplayY:m,detectPointToDisplay:y=>qt(y,p,m),measurePointToDisplay:y=>qt(y,p,m),displayPointToDetect:y=>qt(y,1/Math.max(1e-9,p),1/Math.max(1e-9,m)),measureUsesDisplayLine:!1,measureWidth:t.width,measureHeight:t.height,refineLine:(y,x,b)=>ye(f,t.width,t.height,y,x,b*.5,Math.max(4,b*.2)),measureEdge:(y,x,b,_,k)=>Pe(t.data,t.width,t.height,y,x,b,_,{isThreePlane:!0,threePlaneChannel:l,preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?ge(k==null?void 0:k.fitPoints,1,1):void 0})}}const d=$n(n,!!e.sfrHasGamma);return{sourceMode:i,detectionGray:d,detectionWidth:n.width,detectionHeight:n.height,detectToMeasureX:1,detectToMeasureY:1,detectToDisplayX:1,detectToDisplayY:1,measureToDisplayX:1,measureToDisplayY:1,detectPointToDisplay:f=>f,measurePointToDisplay:f=>f,displayPointToDetect:f=>f,measureUsesDisplayLine:!1,measureWidth:n.width,measureHeight:n.height,refineLine:(f,p,m)=>ye(d,n.width,n.height,f,p,m*.5,Math.max(4,m*.2)),measureEdge:(f,p,m,y,x)=>{const b=$r(n,f,!!e.sfrHasGamma);if(!b)return null;const _=fn(p,f.x,f.y);return Pe(b.data,b.width,b.height,{x:0,y:0,w:b.width,h:b.height},_,m,y,{isThreePlane:!0,threePlaneChannel:l,preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?ge(x==null?void 0:x.fitPoints,1,1,f.x,f.y):void 0})}}}if(i==="unmix-bw"){if(t&&!t.isThreePlane&&e.displaySettings){const f=ja(t,e.displaySettings,e.blackLevel??e.monochromeBlackLevel??void 0);if(f){const p=wl(f),m=o/Math.max(1,t.width),y=a/Math.max(1,t.height);return{sourceMode:i,detectionGray:p,detectionWidth:t.width,detectionHeight:t.height,detectToMeasureX:1,detectToMeasureY:1,detectToDisplayX:m,detectToDisplayY:y,measureToDisplayX:m,measureToDisplayY:y,detectPointToDisplay:x=>qt(x,m,y),measurePointToDisplay:x=>qt(x,m,y),displayPointToDetect:x=>qt(x,1/Math.max(1e-9,m),1/Math.max(1e-9,y)),measureUsesDisplayLine:!1,measureWidth:t.width,measureHeight:t.height,refineLine:(x,b,_)=>ye(p,t.width,t.height,x,b,_*.5,Math.max(4,_*.2)),measureEdge:(x,b,_,k,v)=>Pe(f.data,f.width,f.height,x,b,_,k,{preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?ge(v==null?void 0:v.fitPoints,1,1):void 0})}}}const d=$n(n,!!e.sfrHasGamma,h);return{sourceMode:i,detectionGray:d,detectionWidth:n.width,detectionHeight:n.height,detectToMeasureX:1,detectToMeasureY:1,detectToDisplayX:1,detectToDisplayY:1,measureToDisplayX:1,measureToDisplayY:1,detectPointToDisplay:f=>f,measurePointToDisplay:f=>f,displayPointToDetect:f=>f,measureUsesDisplayLine:!1,measureWidth:n.width,measureHeight:n.height,refineLine:(f,p,m)=>ye(d,n.width,n.height,f,p,m*.5,Math.max(4,m*.2)),measureEdge:(f,p,m,y,x)=>{const b=Jr(n,f,!!e.sfrHasGamma,h);if(!b)return null;const _=fn(p,f.x,f.y);return Pe(b.data,b.width,b.height,{x:0,y:0,w:b.width,h:b.height},_,m,y,{isThreePlane:!0,threePlaneChannel:void 0,preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?ge(x==null?void 0:x.fitPoints,1,1,f.x,f.y):void 0})}}}const c=$n(n,!1);if(t&&!t.isThreePlane&&((g=e.displaySettings)==null?void 0:g.renderMode)==="advanced-zero-dep"&&e.displaySettings.advancedZeroDep){const d=o/Math.max(1,t.width),f=a/Math.max(1,t.height);return{sourceMode:i,detectionGray:c,detectionWidth:n.width,detectionHeight:n.height,detectToMeasureX:t.width/Math.max(1,n.width),detectToMeasureY:t.height/Math.max(1,n.height),detectToDisplayX:1,detectToDisplayY:1,measureToDisplayX:d,measureToDisplayY:f,detectPointToDisplay:p=>p,measurePointToDisplay:p=>qt(p,d,f),displayPointToDetect:p=>p,measureUsesDisplayLine:!1,measureWidth:t.width,measureHeight:t.height,refineLine:(p,m,y)=>ye(c,n.width,n.height,p,m,y*.5,Math.max(4,y*.2)),measureEdge:(p,m,y,x,b)=>{const _=Wa(t,p,e.displaySettings);if(!_||_.width<8||_.height<8)return null;const k=fn(m,p.x,p.y);return Pe(_.data,_.width,_.height,{x:0,y:0,w:_.width,h:_.height},k,y,x,{preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?ge(b==null?void 0:b.fitPoints,t.width/Math.max(1,n.width),t.height/Math.max(1,n.height),p.x,p.y):void 0})}}}if(t&&!t.isThreePlane){const d=o/Math.max(1,t.width),f=a/Math.max(1,t.height);return{sourceMode:i,detectionGray:c,detectionWidth:n.width,detectionHeight:n.height,detectToMeasureX:t.width/Math.max(1,n.width),detectToMeasureY:t.height/Math.max(1,n.height),detectToDisplayX:1,detectToDisplayY:1,measureToDisplayX:d,measureToDisplayY:f,detectPointToDisplay:p=>p,measurePointToDisplay:p=>qt(p,d,f),displayPointToDetect:p=>p,measureUsesDisplayLine:!1,measureWidth:t.width,measureHeight:t.height,refineLine:(p,m,y)=>ye(c,n.width,n.height,p,m,y*.5,Math.max(4,y*.2)),measureEdge:(p,m,y,x,b)=>Pe(t.data,t.width,t.height,p,m,y,x,{blackLevel:e.blackLevel??void 0,preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?ge(b==null?void 0:b.fitPoints,t.width/Math.max(1,n.width),t.height/Math.max(1,n.height)):void 0})}}return{sourceMode:i,detectionGray:c,detectionWidth:n.width,detectionHeight:n.height,detectToMeasureX:1,detectToMeasureY:1,detectToDisplayX:1,detectToDisplayY:1,measureToDisplayX:1,measureToDisplayY:1,detectPointToDisplay:d=>d,measurePointToDisplay:d=>d,displayPointToDetect:d=>d,measureUsesDisplayLine:!1,measureWidth:n.width,measureHeight:n.height,refineLine:(d,f,p)=>ye(c,n.width,n.height,d,f,p*.5,Math.max(4,p*.2)),measureEdge:(d,f,p,m,y)=>{const x=Sl(n,d);if(!x)return null;const b=fn(f,d.x,d.y);return Vo(x.data,x.width,x.height,{x:0,y:0,w:x.width,h:x.height},b,p,m,{blackLevel:e.monochromeBlackLevel??void 0,preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?ge(y==null?void 0:y.fitPoints,1,1,d.x,d.y):void 0})}}}function Xl(n,t,e){var u,h,c,g,d,f,p,m;if(!n||!t)return[];(u=e.onProgress)==null||u.call(e,"Preparing source context...",0);const i=Vl(n,t,e);if(!i)return[];(h=e.onProgress)==null||h.call(e,"Preparing source context...",.08);const r=ar(e.detectionTuning),s=Math.min(1e3,Math.max(1,e.maxRegions??1e3)),o=Math.max(4,e.maxEdges??s*4);(c=e.onProgress)==null||c.call(e,"Detecting candidates...",.12);const a=ll(i.detectionGray,i.detectionWidth,i.detectionHeight,s,e.detectionTuning,(y,x)=>{var b;(b=e.onProgress)==null||b.call(e,y,.12+.08*Math.max(0,Math.min(1,x)))});(g=e.onProgress)==null||g.call(e,"Detecting candidates...",.2);const l=[];for(let y=0;y<a.length;y++){const x=a.length<=0?1:y/a.length;if((d=e.onProgress)==null||d.call(e,`Measuring edges: region ${y+1}/${a.length}`,.2+.72*Math.min(1,x)),l.length>=o)break;const b=a[y],_=b.corners,k=`auto-region-${y+1}`;for(let v=0;v<4&&((f=e.onProgress)==null||f.call(e,`Measuring edges: region ${y+1}/${a.length}, edge ${v+1}/4`,.2+.72*Math.min(1,(y+v/4)/Math.max(1,a.length))),!(l.length>=o));v=v+1){const w=_[v],P=_[(v+1)%4],A=P.x-w.x,F=P.y-w.y,C=Math.hypot(A,F);if(!Number.isFinite(C)||C<24)continue;const M=.125,S={x:w.x+A*M,y:w.y+F*M},I={x:P.x-A*M,y:P.y-F*M},N=Math.hypot(I.x-S.x,I.y-S.y);if(!Number.isFinite(N)||N<12)continue;const E=i.refineLine(S,I,N),X=(E!=null&&E.fitPoints?zr(E.fitPoints):null)||(E==null?void 0:E.line)||{p1:S,p2:I},U=Pl(X,i.detectToMeasureX,i.detectToMeasureY),L=Cl((E==null?void 0:E.fitPoints)??[],i.detectPointToDisplay),B=(L.length>=2?zr(L):null)||kl(U,i.measurePointToDisplay),O=i.measureUsesDisplayLine?B:U,z=O.p2.x-O.p1.x,J=O.p2.y-O.p1.y,V=Math.hypot(z,J);if(!Number.isFinite(V)||V<=1e-6)continue;const q=B.p2.x-B.p1.x,W=B.p2.y-B.p1.y,tt=Math.hypot(q,W);if(!Number.isFinite(tt)||tt<=1e-6)continue;const at=!!e.distortionCurveApplied&&!!e.distortionModel,D=q/tt,Z=W/tt;let H=Z,T=-D;const Q=(B.p1.x+B.p2.x)*.5,G=(B.p1.y+B.p2.y)*.5,Y=i.detectPointToDisplay({x:b.centerX,y:b.centerY}),et=Y.x,nt=Y.y;(Q-et)*H+(G-nt)*T<0&&(H=-H,T=-T);const it=V*.5,ct=Math.max(2,V*r.sampleHalfWidthRatio),dt=Math.max(2,tt*r.sampleHalfWidthRatio),rt=at?{p1:Ut(B.p1,e.distortionModel),p2:Ut(B.p2,e.distortionModel)}:void 0,Bt=rt?Math.max(1,Math.hypot(rt.p2.x-rt.p1.x,rt.p2.y-rt.p1.y)*.5):it,Xt=at?B:U,j=at?dt:ct,ut=_e(Xt,j);if(!ut)continue;const ot=Mt(Vt(ut,2),at?n.width:i.measureWidth,at?n.height:i.measureHeight);if(!ot)continue;const pt=e.distortionCorrected&&e.distortionModel&&i.sourceMode==="rggb-raw"?sr(ot,e.distortionModel,t.width,t.height):null,mt=at?Os(e.distortionSamplingPlane??e.distortionSamplingImage??n,e.distortionModel,rt,Bt,j,!!e.sfrHasGamma,ot,e.distortionBaseImage??n,Xt):i.measureEdge(ot,O,it,j,E);if(!mt||(mt.autoLikeUsed=!0,!Do(mt,e.useDeshading,0)))continue;const ee=e.useNR?-1:12,Dt=Wo([mt],ee,null,e.useDeshading,0,!0);if(!Dt||Dt.mtf50===null||!Oo(Dt.lsfCropped))continue;const _t=rt?Nl(rt,e.distortionModel,Math.max(21,Math.round(tt*.5))):mt.quadraticProjectionUsed?Bo(L,B,Math.max(21,Math.round(tt*.5))):void 0,ce=rt&&_t&&_t.length>=2?Vt(_t,dt+2):null,$t=ce?Rl(ce):_e(B,dt);if(!$t)continue;const Ae=e.distortionCorrected?ot:ce??Vt($t,2);let $={x:Q+H*(dt+12),y:G+T*(dt+12)},Ft=Yr(D,Z);if(_t&&_t.length>=3){const Tt=Math.floor(_t.length/2),xt=_t[Math.max(0,Tt-1)],wt=_t[Math.min(_t.length-1,Tt+1)],ne=_t[Tt],Jt=wt.x-xt.x,Sn=wt.y-xt.y,Oe=Math.hypot(Jt,Sn);if(Oe>1e-6){const Ve=Sn/Oe,Xe=-Jt/Oe;Ft=Yr(Jt/Oe,Sn/Oe);const Ge={x:ne.x-et,y:ne.y-nt},ze=Ge.x*Ve+Ge.y*Xe>=0?1:-1;$={x:ne.x+Ve*ze*(dt+12),y:ne.y+Xe*ze*(dt+12)}}}l.push({id:`${k}-edge-${v+1}`,regionId:k,sourceMode:i.sourceMode,edgeIndex:v,label:Dt.mtf50.toFixed(3),mtf50:Dt.mtf50,angle:Ft,orientation:mt.orientation,edgeData:mt,sourceRect:Ae,rawSourceRect:(i.sourceMode==="rggb-raw"?pt??ot:pt)??void 0,quad:$t,line:B,originalLine:U,curveBaseLine:rt,curvePoints:_t,labelPoint:$,ridgePoints:L,outerSideMeans:b.outerSideMeans,outerSideQuads:b.outerSideQuads,distortionCorrected:e.distortionCorrected??!1})}}return(p=e.onProgress)==null||p.call(e,"Finalizing results...",.98),(m=e.onProgress)==null||m.call(e,"Finalizing results...",1),l}const Gl=n=>!n.blackLevels||n.blackLevels.length<4?null:[Number(n.blackLevels[0])||0,Number(n.blackLevels[1])||0,Number(n.blackLevels[2])||0,Number(n.blackLevels[3])||0],Oi=(n,t)=>{t instanceof ArrayBuffer&&(n.includes(t)||n.push(t))};self.onmessage=async n=>{var o,a;const{id:t,buffer:e,detect:i,options:r}=n.data,s=performance.now();try{const l=performance.now(),u=await Ya(e),h=performance.now()-l;let c=0,g=[];if(i&&!u.isXTrans){const f=u.isThreePlane?"three-plane":"rggb-raw",p=performance.now();g=Xl({width:u.width,height:u.height},u,{...r,sourceMode:f,forceRenderedMeasurement:!1,blackLevel:(r==null?void 0:r.blackLevel)??Gl(u),onProgress:(m,y)=>{self.postMessage({id:t,type:"progress",stage:m,progress:y})}}),c=performance.now()-p}const d=[];Oi(d,e),Oi(d,(o=u.data)==null?void 0:o.buffer),Oi(d,(a=u.floatData)==null?void 0:a.buffer),self.postMessage({id:t,type:"result",success:!0,raw:u,rawFileBuffer:e,measurements:g,timings:{decodeMs:h,detectMs:c,totalMs:performance.now()-s}},d)}catch(l){self.postMessage({id:t,type:"result",success:!1,error:(l==null?void 0:l.message)||String(l)})}};
