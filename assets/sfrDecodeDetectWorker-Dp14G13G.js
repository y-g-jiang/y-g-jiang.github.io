var ra=Object.defineProperty;var sa=(n,t,e)=>t in n?ra(n,t,{enumerable:!0,configurable:!0,writable:!0,value:e}):n[t]=e;var yt=(n,t,e)=>sa(n,typeof t!="symbol"?t+"":t,e);function aa(n){const[[t,e,i],[r,s,o],[a,l,u]]=n,h=t*(s*u-o*l)-e*(r*u-o*a)+i*(r*l-s*a);if(Math.abs(h)<1e-10)return null;const c=1/h;return[[(s*u-o*l)*c,(i*l-e*u)*c,(e*o-i*s)*c],[(o*a-r*u)*c,(t*u-i*a)*c,(i*r-t*o)*c],[(r*l-s*a)*c,(e*a-t*l)*c,(t*s-e*r)*c]]}function oa(n,t){const[[e,i,r],[s,o,a],[l,u,h]]=n,[c,y,f]=t;return[e*c+i*y+r*f,s*c+o*y+a*f,l*c+u*y+h*f]}class gn{constructor(t){yt(this,"size");yt(this,"isPowerOfTwo");yt(this,"_real");yt(this,"_imag");yt(this,"_scratch",null);yt(this,"_rev",null);yt(this,"_m",0);yt(this,"_internalFFT",null);yt(this,"_chirpReal",null);yt(this,"_chirpImag",null);yt(this,"_bReal",null);yt(this,"_bImag",null);yt(this,"_hanning",null);yt(this,"_windowSumSq",0);this.size=t,this.isPowerOfTwo=(t&t-1)===0&&t>0,this._real=new Float32Array(t),this._imag=new Float32Array(t),this.isPowerOfTwo?this.initRadix2():this.initBluestein()}initRadix2(){const t=this.size,e=Math.log2(t);this._rev=new Uint32Array(t);for(let i=0;i<t;i++){let r=0,s=i;for(let o=0;o<e;o++)r=r<<1|s&1,s>>>=1;this._rev[i]=r}}initBluestein(){const t=this.size;this._m=Math.pow(2,Math.ceil(Math.log2(2*t-1))),this._internalFFT=new gn(this._m),this._chirpReal=new Float32Array(t),this._chirpImag=new Float32Array(t);for(let r=0;r<t;r++){const s=-Math.PI*(r*r)/t;this._chirpReal[r]=Math.cos(s),this._chirpImag[r]=Math.sin(s)}const e=new Float32Array(this._m),i=new Float32Array(this._m);for(let r=0;r<t;r++)e[r]=this._chirpReal[r],i[r]=-this._chirpImag[r];for(let r=1;r<t;r++)e[this._m-r]=e[r],i[this._m-r]=i[r];this._internalFFT.transform(e,i),this._bReal=new Float32Array(this._internalFFT._real),this._bImag=new Float32Array(this._internalFFT._imag)}initHanning(){if(this._hanning)return;const t=this.size;this._hanning=new Float32Array(t);let e=0;for(let i=0;i<t;i++){const r=.5*(1-Math.cos(2*Math.PI*i/(t-1)));this._hanning[i]=r,e+=r*r}this._windowSumSq=e}transform(t,e){this.isPowerOfTwo?this.transformRadix2(t,e):this.transformBluestein(t,e)}transformRadix2(t,e){const i=this.size,r=this._rev,s=this._real,o=this._imag;if(t===s)for(let a=0;a<i;a++){const l=r[a];if(a<l){const u=s[a],h=o[a];s[a]=s[l],o[a]=o[l],s[l]=u,o[l]=h}}else for(let a=0;a<i;a++){const l=r[a];s[a]=t[l],o[a]=e?e[l]:0}for(let a=2;a<=i;a*=2){const l=a/2,u=-2*Math.PI/a,h=Math.cos(u),c=Math.sin(u);for(let y=0;y<i;y+=a){let f=1,p=0;for(let d=0;d<l;d++){const g=y+d,x=y+d+l,m=f*s[x]-p*o[x],b=f*o[x]+p*s[x],_=s[g],P=o[g];s[g]=_+m,o[g]=P+b,s[x]=_-m,o[x]=P-b;const M=f*h-p*c,w=f*c+p*h;f=M,p=w}}}}transformBluestein(t,e){const i=this.size,r=this._m,s=this._internalFFT,o=s._real,a=s._imag;o.fill(0),a.fill(0);for(let c=0;c<i;c++){const y=t[c],f=e?e[c]:0,p=this._chirpReal[c],d=this._chirpImag[c];o[c]=y*p-f*d,a[c]=y*d+f*p}s.transformRadix2(o,a);for(let c=0;c<r;c++){const y=s._real[c],f=s._imag[c],p=this._bReal[c],d=this._bImag[c];s._real[c]=y*p-f*d,s._imag[c]=y*d+f*p}const l=s._real,u=s._imag;for(let c=0;c<r;c++)u[c]=-u[c];s.transformRadix2(l,u);const h=1/r;for(let c=0;c<i;c++){const y=s._real[c]*h,f=-s._imag[c]*h,p=this._chirpReal[c],d=this._chirpImag[c];this._real[c]=y*p-f*d,this._imag[c]=y*d+f*p}}calculateSpectrum(t,e,i=!1){const r=this.size;let s=0;for(let h=0;h<r;h++)s+=t[h];const o=s/r;this._scratch||(this._scratch=new Float32Array(r));const a=this._scratch;if(i){this.initHanning();const h=this._hanning;for(let c=0;c<r;c++)a[c]=(t[c]-o)*h[c]}else for(let h=0;h<r;h++)a[h]=t[h]-o;this.transform(a);const l=e.length;let u=1/r;i&&this._windowSumSq>0&&(u=1/this._windowSumSq);for(let h=0;h<l;h++){const c=this._real[h],y=this._imag[h];e[h]+=(c*c+y*y)*u}}calculateSpectrumWindow(t,e,i,r=!1){const s=this.size;let o=0;for(let c=0;c<s;c++)o+=t[e+c];const a=o/s;this._scratch||(this._scratch=new Float32Array(s));const l=this._scratch;if(r){this.initHanning();const c=this._hanning;for(let y=0;y<s;y++)l[y]=(t[e+y]-a)*c[y]}else for(let c=0;c<s;c++)l[c]=t[e+c]-a;this.transform(l);const u=i.length;let h=1/s;r&&this._windowSumSq>0&&(h=1/this._windowSumSq);for(let c=0;c<u;c++){const y=this._real[c],f=this._imag[c];i[c]+=(y*y+f*f)*h}}}const la={"Sony ILCE-7RM5":"0.82 -0.2976 -0.0719 -0.4296 1.2053 0.2532 -0.0429 0.1282 0.5774"};let Ii=null;async function ca(n){return Ii||(Ii=(async()=>{if(typeof window.loadPyodide!="function")throw new Error("Pyodide missing: window.loadPyodide not found.");const t=await window.loadPyodide();return await t.loadPackage("numpy"),t})()),Ii}let Ni=null,Cr=!1;async function $i(){var e;Ni||(Ni=(async()=>{const i=await import("./joraw2-Bb3_vNP4.js");if(typeof i.default!="function")throw new Error("JoRaw2 WASM import failed");const r=new URL("/assets/joraw2-3YkywkGx.wasm",import.meta.url).href;return i.default({locateFile(s,o){return s.endsWith("joraw2.wasm")?r:o+s}})})());const t=(await Ni).LibRaw;if(!t)throw new Error("JoRaw2 class not found");if(!Cr){const i=new t;try{if(typeof i.runtimeInfo!="function")throw new Error("JoRaw2 runtime identity is missing");const r=i.runtimeInfo();if((r==null?void 0:r.wrapper)!=="joraw2"||!String((r==null?void 0:r.librawVersion)||"").startsWith("0.22.2")||!(r!=null&&r.nikonHe)||!(r!=null&&r.nikonHeStar))throw new Error(`Unexpected JoRaw2 runtime: ${JSON.stringify(r)}`);Cr=!0}finally{(e=i.delete)==null||e.call(i)}}return t}var ua=`#!/usr/bin/env python3
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
`,ha=`#!/usr/bin/env python3
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
`,da=`#!/usr/bin/env python3
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
`;const kr=512,fa=la["Sony ILCE-7RM5"].split(/\s+/).map(Number).filter(Number.isFinite);let Ri=null;function pa(n){if(n.byteLength<8)return null;const t=n.getUint16(0,!1);return t===18761?!0:t===19789?!1:null}function hs(n,t,e){const i=Math.min(n.length,t+e);let r="";for(let s=t;s<i;s++){const o=n[s];if(o===0)break;r+=String.fromCharCode(o)}return r.trim()}function Vt(n,t,e,i,r){const s=e===1||e===2||e===7?1:e===3||e===8?2:e===4||e===9?4:0;if(!s)return[];const o=s*i,a=o<=4?r:n.getUint32(r,t);if(a<0||a+o>n.byteLength)return[];const l=[];for(let u=0;u<i;u++){const h=a+u*s;e===1||e===2||e===7?l.push(n.getUint8(h)):e===3?l.push(n.getUint16(h,t)):e===8?l.push(n.getInt16(h,t)):e===4?l.push(n.getUint32(h,t)):e===9&&l.push(n.getInt32(h,t))}return l}function Fr(n,t,e,i,r,s){if(i!==2||r<=0)return"";const o=r<=4?s:t.getUint32(s,e);return o<0||o>=n.length?"":hs(n,o,r)}function ds(n){const t=new Uint8Array(n),e=new DataView(n),i=pa(e);if(i===null||e.getUint16(2,i)!==42)return null;const s=c=>e.getUint16(c,i),o=c=>e.getUint32(c,i),a=[o(4)],l=new Set;let u="",h="";for(;a.length;){const c=a.pop();if(l.has(c)||c<=0||c+2>e.byteLength)continue;l.add(c);const y=s(c);if(c+2+y*12+4>e.byteLength)continue;const f=new Map;for(let _=0;_<y;_++){const P=c+2+_*12,M=s(P),w=s(P+2),C=o(P+4),F=P+8;f.set(M,{type:w,count:C,valueOffset:F})}const p=f.get(271),d=f.get(272);p&&!u&&(u=Fr(t,e,i,p.type,p.count,p.valueOffset)),d&&!h&&(h=Fr(t,e,i,d.type,d.count,d.valueOffset));const g=f.get(330);if(g){const _=Vt(e,i,g.type,g.count,g.valueOffset);for(const P of _)a.push(P)}const x=f.get(259),m=f.get(262);if(x&&m){const _=Vt(e,i,x.type,x.count,x.valueOffset)[0],P=Vt(e,i,m.type,m.count,m.valueOffset)[0];if(_===32766&&P===32803){const M=Vt(e,i,f.get(256).type,f.get(256).count,f.get(256).valueOffset)[0],w=Vt(e,i,f.get(257).type,f.get(257).count,f.get(257).valueOffset)[0],C=Vt(e,i,f.get(258).type,f.get(258).count,f.get(258).valueOffset)[0],F=Vt(e,i,f.get(273).type,f.get(273).count,f.get(273).valueOffset)[0],A=Vt(e,i,f.get(279).type,f.get(279).count,f.get(279).valueOffset)[0],k=f.get(33422)?Vt(e,i,f.get(33422).type,f.get(33422).count,f.get(33422).valueOffset):[0,1,1,2];f.get(29456)&&Vt(e,i,f.get(29456).type,f.get(29456).count,f.get(29456).valueOffset);const v=f.get(50717)?Vt(e,i,f.get(50717).type,f.get(50717).count,f.get(50717).valueOffset)[0]:16383,S=f.get(50719)?Vt(e,i,f.get(50719).type,f.get(50719).count,f.get(50719).valueOffset):[],T=f.get(50720)?Vt(e,i,f.get(50720).type,f.get(50720).count,f.get(50720).valueOffset):[];if(F+kr+16>t.length||F+A>t.length)return null;const L=F+kr,E=hs(t,L,4),V=t[L+8]<<8|t[L+9],U=t[L+10]<<8|t[L+11],N=t[L+12]<<8|t[L+13],D=t[L+14]<<8|t[L+15],B=N>>4&63,z=D>>13,q=D>>10&3,O=U*2,W=A>=4?(t[F]|t[F+1]<<8|t[F+2]<<16|t[F+3]<<24)>>>0:0,j=V===M&&O===w;let K=!1;if(W>=1&&W<=16&&A>=8+W*24){const rt=new Map,R=new Map;let Q=!0;for(let G=0;G<W;G++){const X=F+8+G*24,st=e.getUint32(X+8,!0),$=e.getUint32(X+12,!0),Z=e.getUint32(X+16,!0),lt=e.getUint32(X+20,!0);if(!Z||!lt||st+Z>M||$+lt>w){Q=!1;break}const ut=rt.get($);if(ut!==void 0&&ut!==lt){Q=!1;break}rt.set($,lt),R.set($,(R.get($)||0)+Z)}if(Q){const G=Array.from(rt.keys()).sort((st,$)=>st-$);let X=0;for(const st of G){if(st!==X||R.get(st)!==M){Q=!1;break}X+=rt.get(st)}K=Q&&X===w}}const et=W>=1&&W<=16&&V>0&&O>0&&M%V===0&&w%O===0&&W===M/V*(w/O);if(E!=="A000"&&E!=="0000"||!j&&!et&&!K||B!==16||z!==3||q!==3)return null;const tt=[1024,1024,1024,1024],ot=k.length>=4?k.slice(0,4).map(rt=>rt===0?"R":rt===2?"B":rt===1?"G":"?").join(""):"";return{width:M,height:w,bitsPerSample:C,compression:_,photometric:P,blackLevel:tt,whiteLevel:Number(v||16383),cfaPattern:ot,defaultCropOrigin:S.length>=2?[Number(S[0]),Number(S[1])]:void 0,defaultCropSize:T.length>=2?[Number(T[0]),Number(T[1])]:void 0,make:u||"SONY",model:h||"ILCE-7M5"}}}const b=o(c+2+y*12);b&&a.push(b)}return null}async function ma(n){return Ri||(Ri=(async()=>{const t=await ca();return t.__jtrSonyCrawHqDecoderReady||(await t.FS.mkdirTree("/sony_craw_hq"),await t.FS.writeFile("/sony_craw_hq/llvc3_bitstream_probe.py",ua),await t.FS.writeFile("/sony_craw_hq/llvc3_entropy.py",ha),await t.FS.writeFile("/sony_craw_hq/llvc3_math.py",da),await t.runPythonAsync(`
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
`),t.__jtrSonyCrawHqDecoderReady=!0),t})()),Ri}function ga(n){return ds(n)}function ya(n){return n==="ILCE-7M5"?fa:null}function fs(n,t,e,i){var l;if(n.length!==t.width*t.height)throw new Error(`Sony cRAW HQ decoded size mismatch: got ${n.length}, expected ${t.width*t.height}`);const r=t.model||"ILCE-7M5",s=r.startsWith("Sony ")?r:`Sony ${r}`,o=ya(r),a={...i||{},make:t.make||(i==null?void 0:i.camera_make)||"SONY",model:r,camera_make:t.make||(i==null?void 0:i.camera_make)||"SONY",camera_model:r,UniqueCameraModel:s,sourceFormat:e==="joraw2-wasm"?"Sony cRAW HQ / LLVC3 (JoRaw2 WASM)":"Sony cRAW HQ / LLVC3",sonyCrawHq:{...t,decodeBackend:e},color_desc:t.cfaPattern,black_level_per_channel:t.blackLevel,white_level:t.whiteLevel,color_matrix:o&&o.length===9?o:void 0,idata:{filters:2492765332,colors:3},color_data:{...(i==null?void 0:i.color_data)||{},black:1024,cblack_rawpy_style:t.blackLevel,dng_levels:{...((l=i==null?void 0:i.color_data)==null?void 0:l.dng_levels)||{},dng_cblack:t.blackLevel,dng_whitelevel:t.whiteLevel}}};return{data:n,width:t.width,height:t.height,bayerPattern:t.cfaPattern,blackLevels:t.blackLevel,whiteLevel:t.whiteLevel,metadata:a,isThreePlane:!1,isXTrans:!1}}async function xa(n,t,e){const i=typeof performance<"u"?performance.now():Date.now(),r=await $i(),s=typeof performance<"u"?performance.now():Date.now(),o=new r;try{const a=new Uint8Array(n);await o.open(a,{});const l=typeof performance<"u"?performance.now():Date.now();let u=null;try{u=await o.metadata(!0)}catch(g){console.warn("[Sony cRAW HQ] fast WASM metadata read failed",g)}const h=typeof performance<"u"?performance.now():Date.now(),c=o.getRawImage(),y=typeof performance<"u"?performance.now():Date.now();if(!c||!c.data)throw new Error("Sony cRAW HQ LibRaw WASM returned no raw image");const f=c.data instanceof Uint16Array?c.data:new Uint16Array(c.data.buffer,c.data.byteOffset||0,c.data.byteLength/2),p=typeof performance<"u"?performance.now():Date.now();if(c.width!==t.width||c.height!==t.height)throw new Error(`Sony cRAW HQ LibRaw WASM dimensions mismatch: got ${c.width}x${c.height}, expected ${t.width}x${t.height}`);const d=typeof performance<"u"?performance.now():Date.now();return console.info("[Sony cRAW HQ] fast decode timings",{width:t.width,height:t.height,backend:"joraw2-wasm",wasmReadyMs:Math.round(s-i),openMs:Math.round(l-s),metadataMs:Math.round(h-l),unpackMs:Math.round(y-h),copyMs:Math.round(p-y),totalMs:Math.round(d-i)}),{rawImageData:fs(f,t,"joraw2-wasm",u),info:t}}finally{typeof o.delete=="function"?o.delete():typeof o.close=="function"&&o.close()}}async function ba(n,t,e){const i=typeof performance<"u"?performance.now():Date.now(),r=await ma(),s=typeof performance<"u"?performance.now():Date.now(),o=new Uint8Array(n),a=await fetch(new URL("/assets/sony_llvc3_static_lut4096_padded_u16-FsVBk-IV.bin",import.meta.url));if(!a.ok)throw new Error(`Failed to load Sony LLVC3 sample LUT: HTTP ${a.status}`);const l=new Uint8Array(await a.arrayBuffer()),u=typeof performance<"u"?performance.now():Date.now();r.globals.set("jtr_sony_arw_bytes",o),r.globals.set("jtr_sony_lut_bytes",l);const h=await r.runPythonAsync("jtr_decode_sony_craw_hq(jtr_sony_arw_bytes.to_py(), jtr_sony_lut_bytes.to_py())"),c=typeof performance<"u"?performance.now():Date.now(),y=h.toJs();typeof h.destroy=="function"&&h.destroy(),r.globals.delete("jtr_sony_arw_bytes"),r.globals.delete("jtr_sony_lut_bytes");const f=new Uint8Array(y.byteLength);f.set(y);const p=new Uint16Array(f.buffer),d=typeof performance<"u"?performance.now():Date.now(),g=typeof performance<"u"?performance.now():Date.now();return console.info("[Sony cRAW HQ] decode timings",{width:t.width,height:t.height,backend:"pyodide",pyodideReadyMs:Math.round(s-i),lutLoadMs:Math.round(u-s),llvc3DecodeMs:Math.round(c-u),copyMs:Math.round(d-c),totalMs:Math.round(g-i)}),{rawImageData:fs(p,t,"pyodide"),info:t}}async function _a(n,t){const e=ds(n);if(!e)return null;try{return await xa(n,e,t)}catch(i){return console.warn("[Sony cRAW HQ] fast WASM decode failed; falling back to Pyodide",i),ba(n,e)}}async function wa(n,t){return _a(n,t)}const Ma=["RGGB","BGGR","GRBG","GBRG"],Sa=new Set(Ma);function va(n){if(!n||typeof n!="object"||Array.isArray(n)||ArrayBuffer.isView(n))return n;const t=n;return t.value??t.values??t.description??n}function Ke(n){const t=va(n);let e="";if(typeof t=="string")e=t.toUpperCase().replace(/[^RGB012]/g,"");else if(typeof t=="number")e=String(t);else if(Array.isArray(t)||ArrayBuffer.isView(t))e=Array.from(t).map(String).join("");else return null;return/^[012]{4}$/.test(e)&&(e=e.replace(/0/g,"R").replace(/1/g,"G").replace(/2/g,"B")),Sa.has(e)?e:null}function Ar(n){var t;return n?[n.cfa_pattern,n.cfaPattern,n.BayerPattern,n.CFAPattern2,n.CFAPattern,(t=n.idata)==null?void 0:t.cfa_pattern]:[]}function Pa(n,t,e){const i=Ke(n.bayerPattern),r=n.bayerPatternSource==="manual";if(i&&!r)return{pattern:i,source:n.bayerPatternSource||"decoder"};if(!r)for(const o of Ar(n.metadata)){const a=Ke(o);if(a)return{pattern:a,source:"libraw-metadata"}}for(const o of Ar(t)){const a=Ke(o);if(a)return{pattern:a,source:"metadata"}}const s=Ke(e)||(r?i:null);return s?{pattern:s,source:"manual"}:{pattern:null,source:null}}function Ca(n,t){return n.bayerPattern=t.pattern||"",n.bayerPatternSource=t.source||void 0,t.pattern&&t.source!=="manual"&&(n.metadata={...n.metadata||{},cfa_pattern:t.pattern,cfaPattern:t.pattern}),t}function ps(n,t){return(t&1)<<1|n&1}function ka(n,t,e){const i=Ke(n);if(!i)throw new Error("Bayer CFA pattern is unresolved.");return i[ps(t,e)]}const le=n=>{const t=Number(n);return Number.isFinite(t)?Math.max(0,t):0};function Qe(n){if(!n||typeof n.length!="number")return null;const t=Array.from(n);return t.length<4?null:[le(t[0]),le(t[1]),le(t[2]),le(t[3])]}function Tr(n){if(!n||n.source!=="libraw")return null;const t=Qe(n.channelOffsets),e=Qe(n.channelLevels),i=Qe(n.siteColorIndices),r=Qe(n.siteBaseLevels),s=Qe(n.siteLevels);if(!t||!e||!i||!r||!s)return null;const o=Math.max(0,Math.floor(le(n.repeatRows))),a=Math.max(0,Math.floor(le(n.repeatCols))),l=o*a,u=n.repeatValues&&typeof n.repeatValues.length=="number"?Array.from(n.repeatValues,le):[],h=l>0&&l<=4098&&u.length>=l;return{source:"libraw",common:le(n.common),channelOffsets:t,channelLevels:e,siteColorIndices:i,siteBaseLevels:r,siteLevels:s,repeatRows:h?o:0,repeatCols:h?a:0,repeatOriginY:Math.max(0,Math.floor(le(n.repeatOriginY))),repeatOriginX:Math.max(0,Math.floor(le(n.repeatOriginX))),repeatValues:h?u.slice(0,l):[]}}var ui=typeof self<"u"?self:global;const Nn=typeof navigator<"u",Fa=Nn&&typeof HTMLImageElement>"u",ti=!(typeof global>"u"||typeof process>"u"||!process.versions||!process.versions.node),hi=ui.Buffer,Wn=ui.BigInt,di=!!hi,Aa=n=>n;function ei(n,t=Aa){if(ti)try{return typeof require=="function"?Promise.resolve(t(require(n))):import(n).then(t)}catch{console.warn(`Couldn't load ${n}`)}}let Ji=ui.fetch;const Ta=n=>Ji=n;if(!ui.fetch){const n=ei("http",i=>i),t=ei("https",i=>i),e=(i,{headers:r}={})=>new Promise(async(s,o)=>{let{port:a,hostname:l,pathname:u,protocol:h,search:c}=new URL(i);const y={method:"GET",hostname:l,path:encodeURI(u)+c,headers:r};a!==""&&(y.port=Number(a));const f=(h==="https:"?await t:await n).request(y,p=>{if(p.statusCode===301||p.statusCode===302){let d=new URL(p.headers.location,i).toString();return e(d,{headers:r}).then(s).catch(o)}s({status:p.statusCode,arrayBuffer:()=>new Promise(d=>{let g=[];p.on("data",x=>g.push(x)),p.on("end",()=>d(Buffer.concat(g)))})})});f.on("error",o),f.end()});Ta(e)}function J(n,t,e){return t in n?Object.defineProperty(n,t,{value:e,enumerable:!0,configurable:!0,writable:!0}):n[t]=e,n}const ni=n=>ms(n)?void 0:n,Ia=n=>n!==void 0;function ms(n){return n===void 0||(n instanceof Map?n.size===0:Object.values(n).filter(Ia).length===0)}function St(n){let t=new Error(n);throw delete t.stack,t}function $e(n){return(n=function(t){for(;t.endsWith("\0");)t=t.slice(0,-1);return t}(n).trim())===""?void 0:n}function Xi(n){let t=function(e){let i=0;return e.ifd0.enabled&&(i+=1024),e.exif.enabled&&(i+=2048),e.makerNote&&(i+=2048),e.userComment&&(i+=1024),e.gps.enabled&&(i+=512),e.interop.enabled&&(i+=100),e.ifd1.enabled&&(i+=1024),i+2048}(n);return n.jfif.enabled&&(t+=50),n.xmp.enabled&&(t+=2e4),n.iptc.enabled&&(t+=14e3),n.icc.enabled&&(t+=6e3),t}const Yi=n=>String.fromCharCode.apply(null,n),Ir=typeof TextDecoder<"u"?new TextDecoder("utf-8"):void 0;function gs(n){return Ir?Ir.decode(n):di?Buffer.from(n).toString("utf8"):decodeURIComponent(escape(Yi(n)))}class Dt{static from(t,e){return t instanceof this&&t.le===e?t:new Dt(t,void 0,void 0,e)}constructor(t,e=0,i,r){if(typeof r=="boolean"&&(this.le=r),Array.isArray(t)&&(t=new Uint8Array(t)),t===0)this.byteOffset=0,this.byteLength=0;else if(t instanceof ArrayBuffer){i===void 0&&(i=t.byteLength-e);let s=new DataView(t,e,i);this._swapDataView(s)}else if(t instanceof Uint8Array||t instanceof DataView||t instanceof Dt){i===void 0&&(i=t.byteLength-e),(e+=t.byteOffset)+i>t.byteOffset+t.byteLength&&St("Creating view outside of available memory in ArrayBuffer");let s=new DataView(t.buffer,e,i);this._swapDataView(s)}else if(typeof t=="number"){let s=new DataView(new ArrayBuffer(t));this._swapDataView(s)}else St("Invalid input argument for BufferView: "+t)}_swapArrayBuffer(t){this._swapDataView(new DataView(t))}_swapBuffer(t){this._swapDataView(new DataView(t.buffer,t.byteOffset,t.byteLength))}_swapDataView(t){this.dataView=t,this.buffer=t.buffer,this.byteOffset=t.byteOffset,this.byteLength=t.byteLength}_lengthToEnd(t){return this.byteLength-t}set(t,e,i=Dt){return t instanceof DataView||t instanceof Dt?t=new Uint8Array(t.buffer,t.byteOffset,t.byteLength):t instanceof ArrayBuffer&&(t=new Uint8Array(t)),t instanceof Uint8Array||St("BufferView.set(): Invalid data argument."),this.toUint8().set(t,e),new i(this,e,t.byteLength)}subarray(t,e){return e=e||this._lengthToEnd(t),new Dt(this,t,e)}toUint8(){return new Uint8Array(this.buffer,this.byteOffset,this.byteLength)}getUint8Array(t,e){return new Uint8Array(this.buffer,this.byteOffset+t,e)}getString(t=0,e=this.byteLength){return gs(this.getUint8Array(t,e))}getLatin1String(t=0,e=this.byteLength){let i=this.getUint8Array(t,e);return Yi(i)}getUnicodeString(t=0,e=this.byteLength){const i=[];for(let r=0;r<e&&t+r<this.byteLength;r+=2)i.push(this.getUint16(t+r));return Yi(i)}getInt8(t){return this.dataView.getInt8(t)}getUint8(t){return this.dataView.getUint8(t)}getInt16(t,e=this.le){return this.dataView.getInt16(t,e)}getInt32(t,e=this.le){return this.dataView.getInt32(t,e)}getUint16(t,e=this.le){return this.dataView.getUint16(t,e)}getUint32(t,e=this.le){return this.dataView.getUint32(t,e)}getFloat32(t,e=this.le){return this.dataView.getFloat32(t,e)}getFloat64(t,e=this.le){return this.dataView.getFloat64(t,e)}getFloat(t,e=this.le){return this.dataView.getFloat32(t,e)}getDouble(t,e=this.le){return this.dataView.getFloat64(t,e)}getUintBytes(t,e,i){switch(e){case 1:return this.getUint8(t,i);case 2:return this.getUint16(t,i);case 4:return this.getUint32(t,i);case 8:return this.getUint64&&this.getUint64(t,i)}}getUint(t,e,i){switch(e){case 8:return this.getUint8(t,i);case 16:return this.getUint16(t,i);case 32:return this.getUint32(t,i);case 64:return this.getUint64&&this.getUint64(t,i)}}toString(t){return this.dataView.toString(t,this.constructor.name)}ensureChunk(){}}function Wi(n,t){St(`${n} '${t}' was not loaded, try using full build of exifr.`)}class Zi extends Map{constructor(t){super(),this.kind=t}get(t,e){return this.has(t)||Wi(this.kind,t),e&&(t in e||function(i,r){St(`Unknown ${i} '${r}'.`)}(this.kind,t),e[t].enabled||Wi(this.kind,t)),super.get(t)}keyList(){return Array.from(this.keys())}}var ue=new Zi("file parser"),Mt=new Zi("segment parser"),de=new Zi("file reader");function Na(n,t){return typeof n=="string"?Nr(n,t):Nn&&!Fa&&n instanceof HTMLImageElement?Nr(n.src,t):n instanceof Uint8Array||n instanceof ArrayBuffer||n instanceof DataView?new Dt(n):Nn&&n instanceof Blob?ji(n,t,"blob",sn):void St("Invalid input argument")}function Nr(n,t){return(e=n).startsWith("data:")||e.length>1e4?Hi(n,t,"base64"):ti&&n.includes("://")?ji(n,t,"url",rn):ti?Hi(n,t,"fs"):Nn?ji(n,t,"url",rn):void St("Invalid input argument");var e}async function ji(n,t,e,i){return de.has(e)?Hi(n,t,e):i?async function(r,s){let o=await s(r);return new Dt(o)}(n,i):void St(`Parser ${e} is not loaded`)}async function Hi(n,t,e){let i=new(de.get(e))(n,t);return await i.read(),i}const rn=n=>Ji(n).then(t=>t.arrayBuffer()),sn=n=>new Promise((t,e)=>{let i=new FileReader;i.onloadend=()=>t(i.result||new ArrayBuffer),i.onerror=e,i.readAsArrayBuffer(n)});class Ra extends Map{get tagKeys(){return this.allKeys||(this.allKeys=Array.from(this.keys())),this.allKeys}get tagValues(){return this.allValues||(this.allValues=Array.from(this.values())),this.allValues}}function wt(n,t,e){let i=new Ra;for(let[r,s]of e)i.set(r,s);if(Array.isArray(t))for(let r of t)n.set(r,i);else n.set(t,i);return i}function an(n,t,e){let i,r=n.get(t);for(i of e)r.set(i[0],i[1])}const Ct=new Map,Jt=new Map,Ie=new Map,Pe=["chunked","firstChunkSize","firstChunkSizeNode","firstChunkSizeBrowser","chunkSize","chunkLimit"],yn=["jfif","xmp","icc","iptc","ihdr"],on=["tiff",...yn],bt=["ifd0","ifd1","exif","gps","interop"],Ce=[...on,...bt],ke=["makerNote","userComment"],xn=["translateKeys","translateValues","reviveValues","multiSegment"],Fe=[...xn,"sanitize","mergeOutput","silentErrors"];class ys{get translate(){return this.translateKeys||this.translateValues||this.reviveValues}}class kn extends ys{get needed(){return this.enabled||this.deps.size>0}constructor(t,e,i,r){if(super(),J(this,"enabled",!1),J(this,"skip",new Set),J(this,"pick",new Set),J(this,"deps",new Set),J(this,"translateKeys",!1),J(this,"translateValues",!1),J(this,"reviveValues",!1),this.key=t,this.enabled=e,this.parse=this.enabled,this.applyInheritables(r),this.canBeFiltered=bt.includes(t),this.canBeFiltered&&(this.dict=Ct.get(t)),i!==void 0)if(Array.isArray(i))this.parse=this.enabled=!0,this.canBeFiltered&&i.length>0&&this.translateTagSet(i,this.pick);else if(typeof i=="object"){if(this.enabled=!0,this.parse=i.parse!==!1,this.canBeFiltered){let{pick:s,skip:o}=i;s&&s.length>0&&this.translateTagSet(s,this.pick),o&&o.length>0&&this.translateTagSet(o,this.skip)}this.applyInheritables(i)}else i===!0||i===!1?this.parse=this.enabled=i:St(`Invalid options argument: ${i}`)}applyInheritables(t){let e,i;for(e of xn)i=t[e],i!==void 0&&(this[e]=i)}translateTagSet(t,e){if(this.dict){let i,r,{tagKeys:s,tagValues:o}=this.dict;for(i of t)typeof i=="string"?(r=o.indexOf(i),r===-1&&(r=s.indexOf(Number(i))),r!==-1&&e.add(Number(s[r]))):e.add(i)}else for(let i of t)e.add(i)}finalizeFilters(){!this.enabled&&this.deps.size>0?(this.enabled=!0,ii(this.pick,this.deps)):this.enabled&&this.pick.size>0&&ii(this.pick,this.deps)}}var It={jfif:!1,tiff:!0,xmp:!1,icc:!1,iptc:!1,ifd0:!0,ifd1:!1,exif:!0,gps:!0,interop:!1,ihdr:void 0,makerNote:!1,userComment:!1,multiSegment:!1,skip:[],pick:[],translateKeys:!0,translateValues:!0,reviveValues:!0,sanitize:!0,mergeOutput:!0,silentErrors:!0,chunked:!0,firstChunkSize:void 0,firstChunkSizeNode:512,firstChunkSizeBrowser:65536,chunkSize:65536,chunkLimit:5},Rr=new Map;class ln extends ys{static useCached(t){let e=Rr.get(t);return e!==void 0||(e=new this(t),Rr.set(t,e)),e}constructor(t){super(),t===!0?this.setupFromTrue():t===void 0?this.setupFromUndefined():Array.isArray(t)?this.setupFromArray(t):typeof t=="object"?this.setupFromObject(t):St(`Invalid options argument ${t}`),this.firstChunkSize===void 0&&(this.firstChunkSize=Nn?this.firstChunkSizeBrowser:this.firstChunkSizeNode),this.mergeOutput&&(this.ifd1.enabled=!1),this.filterNestedSegmentTags(),this.traverseTiffDependencyTree(),this.checkLoadedPlugins()}setupFromUndefined(){let t;for(t of Pe)this[t]=It[t];for(t of Fe)this[t]=It[t];for(t of ke)this[t]=It[t];for(t of Ce)this[t]=new kn(t,It[t],void 0,this)}setupFromTrue(){let t;for(t of Pe)this[t]=It[t];for(t of Fe)this[t]=It[t];for(t of ke)this[t]=!0;for(t of Ce)this[t]=new kn(t,!0,void 0,this)}setupFromArray(t){let e;for(e of Pe)this[e]=It[e];for(e of Fe)this[e]=It[e];for(e of ke)this[e]=It[e];for(e of Ce)this[e]=new kn(e,!1,void 0,this);this.setupGlobalFilters(t,void 0,bt)}setupFromObject(t){let e;for(e of(bt.ifd0=bt.ifd0||bt.image,bt.ifd1=bt.ifd1||bt.thumbnail,Object.assign(this,t),Pe))this[e]=Li(t[e],It[e]);for(e of Fe)this[e]=Li(t[e],It[e]);for(e of ke)this[e]=Li(t[e],It[e]);for(e of on)this[e]=new kn(e,It[e],t[e],this);for(e of bt)this[e]=new kn(e,It[e],t[e],this.tiff);this.setupGlobalFilters(t.pick,t.skip,bt,Ce),t.tiff===!0?this.batchEnableWithBool(bt,!0):t.tiff===!1?this.batchEnableWithUserValue(bt,t):Array.isArray(t.tiff)?this.setupGlobalFilters(t.tiff,void 0,bt):typeof t.tiff=="object"&&this.setupGlobalFilters(t.tiff.pick,t.tiff.skip,bt)}batchEnableWithBool(t,e){for(let i of t)this[i].enabled=e}batchEnableWithUserValue(t,e){for(let i of t){let r=e[i];this[i].enabled=r!==!1&&r!==void 0}}setupGlobalFilters(t,e,i,r=i){if(t&&t.length){for(let o of r)this[o].enabled=!1;let s=Lr(t,i);for(let[o,a]of s)ii(this[o].pick,a),this[o].enabled=!0}else if(e&&e.length){let s=Lr(e,i);for(let[o,a]of s)ii(this[o].skip,a)}}filterNestedSegmentTags(){let{ifd0:t,exif:e,xmp:i,iptc:r,icc:s}=this;this.makerNote?e.deps.add(37500):e.skip.add(37500),this.userComment?e.deps.add(37510):e.skip.add(37510),i.enabled||t.skip.add(700),r.enabled||t.skip.add(33723),s.enabled||t.skip.add(34675)}traverseTiffDependencyTree(){let{ifd0:t,exif:e,gps:i,interop:r}=this;r.needed&&(e.deps.add(40965),t.deps.add(40965)),e.needed&&t.deps.add(34665),i.needed&&t.deps.add(34853),this.tiff.enabled=bt.some(s=>this[s].enabled===!0)||this.makerNote||this.userComment;for(let s of bt)this[s].finalizeFilters()}get onlyTiff(){return!yn.map(t=>this[t].enabled).some(t=>t===!0)&&this.tiff.enabled}checkLoadedPlugins(){for(let t of on)this[t].enabled&&!Mt.has(t)&&Wi("segment parser",t)}}function Lr(n,t){let e,i,r,s,o=[];for(r of t){for(s of(e=Ct.get(r),i=[],e))(n.includes(s[0])||n.includes(s[1]))&&i.push(s[0]);i.length&&o.push([r,i])}return o}function Li(n,t){return n!==void 0?n:t!==void 0?t:void 0}function ii(n,t){for(let e of t)n.add(e)}J(ln,"default",It);class Ne{constructor(t){J(this,"parsers",{}),J(this,"output",{}),J(this,"errors",[]),J(this,"pushToErrors",e=>this.errors.push(e)),this.options=ln.useCached(t)}async read(t){this.file=await Na(t,this.options)}setup(){if(this.fileParser)return;let{file:t}=this,e=t.getUint16(0);for(let[i,r]of ue)if(r.canHandle(t,e))return this.fileParser=new r(this.options,this.file,this.parsers),t[i]=!0;this.file.close&&this.file.close(),St("Unknown file format")}async parse(){let{output:t,errors:e}=this;return this.setup(),this.options.silentErrors?(await this.executeParsers().catch(this.pushToErrors),e.push(...this.fileParser.errors)):await this.executeParsers(),this.file.close&&this.file.close(),this.options.silentErrors&&e.length>0&&(t.errors=e),ni(t)}async executeParsers(){let{output:t}=this;await this.fileParser.parse();let e=Object.values(this.parsers).map(async i=>{let r=await i.parse();i.assignToOutput(t,r)});this.options.silentErrors&&(e=e.map(i=>i.catch(this.pushToErrors))),await Promise.all(e)}async extractThumbnail(){this.setup();let{options:t,file:e}=this,i=Mt.get("tiff",t);var r;if(e.tiff?r={start:0,type:"tiff"}:e.jpeg&&(r=await this.fileParser.getOrFindSegment("tiff")),r===void 0)return;let s=await this.fileParser.ensureSegmentChunk(r),o=this.parsers.tiff=new i(s,t,e),a=await o.extractThumbnail();return e.close&&e.close(),a}}async function fi(n,t){let e=new Ne(t);return await e.read(n),e.parse()}var La=Object.freeze({__proto__:null,parse:fi,Exifr:Ne,fileParsers:ue,segmentParsers:Mt,fileReaders:de,tagKeys:Ct,tagValues:Jt,tagRevivers:Ie,createDictionary:wt,extendDictionary:an,fetchUrlAsArrayBuffer:rn,readBlobAsArrayBuffer:sn,chunkedProps:Pe,otherSegments:yn,segments:on,tiffBlocks:bt,segmentsAndBlocks:Ce,tiffExtractables:ke,inheritables:xn,allFormatters:Fe,Options:ln});class pi{constructor(t,e,i){J(this,"errors",[]),J(this,"ensureSegmentChunk",async r=>{let s=r.start,o=r.size||65536;if(this.file.chunked)if(this.file.available(s,o))r.chunk=this.file.subarray(s,o);else try{r.chunk=await this.file.readChunk(s,o)}catch(a){St(`Couldn't read segment: ${JSON.stringify(r)}. ${a.message}`)}else this.file.byteLength>s+o?r.chunk=this.file.subarray(s,o):r.size===void 0?r.chunk=this.file.subarray(s):St("Segment unreachable: "+JSON.stringify(r));return r.chunk}),this.extendOptions&&this.extendOptions(t),this.options=t,this.file=e,this.parsers=i}injectSegment(t,e){this.options[t].enabled&&this.createParser(t,e)}createParser(t,e){let i=new(Mt.get(t))(e,this.options,this.file);return this.parsers[t]=i}createParsers(t){for(let e of t){let{type:i,chunk:r}=e,s=this.options[i];if(s&&s.enabled){let o=this.parsers[i];o&&o.append||o||this.createParser(i,r)}}}async readSegments(t){let e=t.map(this.ensureSegmentChunk);await Promise.all(e)}}class $t{static findPosition(t,e){let i=t.getUint16(e+2)+2,r=typeof this.headerLength=="function"?this.headerLength(t,e,i):this.headerLength,s=e+r,o=i-r;return{offset:e,length:i,headerLength:r,start:s,size:o,end:s+o}}static parse(t,e={}){return new this(t,new ln({[this.type]:e}),t).parse()}normalizeInput(t){return t instanceof Dt?t:new Dt(t)}constructor(t,e={},i){J(this,"errors",[]),J(this,"raw",new Map),J(this,"handleError",r=>{if(!this.options.silentErrors)throw r;this.errors.push(r.message)}),this.chunk=this.normalizeInput(t),this.file=i,this.type=this.constructor.type,this.globalOptions=this.options=e,this.localOptions=e[this.type],this.canTranslate=this.localOptions&&this.localOptions.translate}translate(){this.canTranslate&&(this.translated=this.translateBlock(this.raw,this.type))}get output(){return this.translated?this.translated:this.raw?Object.fromEntries(this.raw):void 0}translateBlock(t,e){let i=Ie.get(e),r=Jt.get(e),s=Ct.get(e),o=this.options[e],a=o.reviveValues&&!!i,l=o.translateValues&&!!r,u=o.translateKeys&&!!s,h={};for(let[c,y]of t)a&&i.has(c)?y=i.get(c)(y):l&&r.has(c)&&(y=this.translateValue(y,r.get(c))),u&&s.has(c)&&(c=s.get(c)||c),h[c]=y;return h}translateValue(t,e){return e[t]||e.DEFAULT||t}assignToOutput(t,e){this.assignObjectToOutput(t,this.constructor.type,e)}assignObjectToOutput(t,e,i){if(this.globalOptions.mergeOutput)return Object.assign(t,i);t[e]?Object.assign(t[e],i):t[e]=i}}J($t,"headerLength",4),J($t,"type",void 0),J($t,"multiSegment",!1),J($t,"canHandle",()=>!1);function Ea(n){return n===192||n===194||n===196||n===219||n===221||n===218||n===254}function Ua(n){return n>=224&&n<=239}function Da(n,t,e){for(let[i,r]of Mt)if(r.canHandle(n,t,e))return i}class Er extends pi{constructor(...t){super(...t),J(this,"appSegments",[]),J(this,"jpegSegments",[]),J(this,"unknownSegments",[])}static canHandle(t,e){return e===65496}async parse(){await this.findAppSegments(),await this.readSegments(this.appSegments),this.mergeMultiSegments(),this.createParsers(this.mergedAppSegments||this.appSegments)}setupSegmentFinderArgs(t){t===!0?(this.findAll=!0,this.wanted=new Set(Mt.keyList())):(t=t===void 0?Mt.keyList().filter(e=>this.options[e].enabled):t.filter(e=>this.options[e].enabled&&Mt.has(e)),this.findAll=!1,this.remaining=new Set(t),this.wanted=new Set(t)),this.unfinishedMultiSegment=!1}async findAppSegments(t=0,e){this.setupSegmentFinderArgs(e);let{file:i,findAll:r,wanted:s,remaining:o}=this;if(!r&&this.file.chunked&&(r=Array.from(s).some(a=>{let l=Mt.get(a),u=this.options[a];return l.multiSegment&&u.multiSegment}),r&&await this.file.readWhole()),t=this.findAppSegmentsInRange(t,i.byteLength),!this.options.onlyTiff&&i.chunked){let a=!1;for(;o.size>0&&!a&&(i.canReadNextChunk||this.unfinishedMultiSegment);){let{nextChunkOffset:l}=i,u=this.appSegments.some(h=>!this.file.available(h.offset||h.start,h.length||h.size));if(a=t>l&&!u?!await i.readNextChunk(t):!await i.readNextChunk(l),(t=this.findAppSegmentsInRange(t,i.byteLength))===void 0)return}}}findAppSegmentsInRange(t,e){e-=2;let i,r,s,o,a,l,{file:u,findAll:h,wanted:c,remaining:y,options:f}=this;for(;t<e;t++)if(u.getUint8(t)===255){if(i=u.getUint8(t+1),Ua(i)){if(r=u.getUint16(t+2),s=Da(u,t,r),s&&c.has(s)&&(o=Mt.get(s),a=o.findPosition(u,t),l=f[s],a.type=s,this.appSegments.push(a),!h&&(o.multiSegment&&l.multiSegment?(this.unfinishedMultiSegment=a.chunkNumber<a.chunkCount,this.unfinishedMultiSegment||y.delete(s)):y.delete(s),y.size===0)))break;f.recordUnknownSegments&&(a=$t.findPosition(u,t),a.marker=i,this.unknownSegments.push(a)),t+=r+1}else if(Ea(i)){if(r=u.getUint16(t+2),i===218&&f.stopAfterSos!==!1)return;f.recordJpegSegments&&this.jpegSegments.push({offset:t,length:r,marker:i}),t+=r+1}}return t}mergeMultiSegments(){if(!this.appSegments.some(e=>e.multiSegment))return;let t=function(e,i){let r,s,o,a=new Map;for(let l=0;l<e.length;l++)r=e[l],s=r[i],a.has(s)?o=a.get(s):a.set(s,o=[]),o.push(r);return Array.from(a)}(this.appSegments,"type");this.mergedAppSegments=t.map(([e,i])=>{let r=Mt.get(e,this.options);return r.handleMultiSegments?{type:e,chunk:r.handleMultiSegments(i)}:i[0]})}getSegment(t){return this.appSegments.find(e=>e.type===t)}async getOrFindSegment(t){let e=this.getSegment(t);return e===void 0&&(await this.findAppSegments(0,[t]),e=this.getSegment(t)),e}}J(Er,"type","jpeg"),ue.set("jpeg",Er);const Ba=[void 0,1,1,2,4,8,1,1,2,4,8,4,8,4];class Oa extends $t{parseHeader(){var t=this.chunk.getUint16();t===18761?this.le=!0:t===19789&&(this.le=!1),this.chunk.le=this.le,this.headerParsed=!0}parseTags(t,e,i=new Map){let{pick:r,skip:s}=this.options[e];r=new Set(r);let o=r.size>0,a=s.size===0,l=this.chunk.getUint16(t);t+=2;for(let u=0;u<l;u++){let h=this.chunk.getUint16(t);if(o){if(r.has(h)&&(i.set(h,this.parseTag(t,h,e)),r.delete(h),r.size===0))break}else!a&&s.has(h)||i.set(h,this.parseTag(t,h,e));t+=12}return i}parseTag(t,e,i){let{chunk:r}=this,s=r.getUint16(t+2),o=r.getUint32(t+4),a=Ba[s];if(a*o<=4?t+=8:t=r.getUint32(t+8),(s<1||s>13)&&St(`Invalid TIFF value type. block: ${i.toUpperCase()}, tag: ${e.toString(16)}, type: ${s}, offset ${t}`),t>r.byteLength&&St(`Invalid TIFF value offset. block: ${i.toUpperCase()}, tag: ${e.toString(16)}, type: ${s}, offset ${t} is outside of chunk size ${r.byteLength}`),s===1)return r.getUint8Array(t,o);if(s===2)return $e(r.getString(t,o));if(s===7)return r.getUint8Array(t,o);if(o===1)return this.parseTagValue(s,t);{let l=new(function(h){switch(h){case 1:return Uint8Array;case 3:return Uint16Array;case 4:return Uint32Array;case 5:return Array;case 6:return Int8Array;case 8:return Int16Array;case 9:return Int32Array;case 10:return Array;case 11:return Float32Array;case 12:return Float64Array;default:return Array}}(s))(o),u=a;for(let h=0;h<o;h++)l[h]=this.parseTagValue(s,t),t+=u;return l}}parseTagValue(t,e){let{chunk:i}=this;switch(t){case 1:return i.getUint8(e);case 3:return i.getUint16(e);case 4:return i.getUint32(e);case 5:return i.getUint32(e)/i.getUint32(e+4);case 6:return i.getInt8(e);case 8:return i.getInt16(e);case 9:return i.getInt32(e);case 10:return i.getInt32(e)/i.getInt32(e+4);case 11:return i.getFloat(e);case 12:return i.getDouble(e);case 13:return i.getUint32(e);default:St(`Invalid tiff type ${t}`)}}}class Ei extends Oa{static canHandle(t,e){return t.getUint8(e+1)===225&&t.getUint32(e+4)===1165519206&&t.getUint16(e+8)===0}async parse(){this.parseHeader();let{options:t}=this;return t.ifd0.enabled&&await this.parseIfd0Block(),t.exif.enabled&&await this.safeParse("parseExifBlock"),t.gps.enabled&&await this.safeParse("parseGpsBlock"),t.interop.enabled&&await this.safeParse("parseInteropBlock"),t.ifd1.enabled&&await this.safeParse("parseThumbnailBlock"),this.createOutput()}safeParse(t){let e=this[t]();return e.catch!==void 0&&(e=e.catch(this.handleError)),e}findIfd0Offset(){this.ifd0Offset===void 0&&(this.ifd0Offset=this.chunk.getUint32(4))}findIfd1Offset(){if(this.ifd1Offset===void 0){this.findIfd0Offset();let t=this.chunk.getUint16(this.ifd0Offset),e=this.ifd0Offset+2+12*t;this.ifd1Offset=this.chunk.getUint32(e)}}parseBlock(t,e){let i=new Map;return this[e]=i,this.parseTags(t,e,i),i}async parseIfd0Block(){if(this.ifd0)return;let{file:t}=this;this.findIfd0Offset(),this.ifd0Offset<8&&St("Malformed EXIF data"),!t.chunked&&this.ifd0Offset>t.byteLength&&St(`IFD0 offset points to outside of file.
this.ifd0Offset: ${this.ifd0Offset}, file.byteLength: ${t.byteLength}`),t.tiff&&await t.ensureChunk(this.ifd0Offset,Xi(this.options));let e=this.parseBlock(this.ifd0Offset,"ifd0");return e.size!==0?(this.exifOffset=e.get(34665),this.interopOffset=e.get(40965),this.gpsOffset=e.get(34853),this.xmp=e.get(700),this.iptc=e.get(33723),this.icc=e.get(34675),this.options.sanitize&&(e.delete(34665),e.delete(40965),e.delete(34853),e.delete(700),e.delete(33723),e.delete(34675)),e):void 0}async parseExifBlock(){if(this.exif||(this.ifd0||await this.parseIfd0Block(),this.exifOffset===void 0))return;this.file.tiff&&await this.file.ensureChunk(this.exifOffset,Xi(this.options));let t=this.parseBlock(this.exifOffset,"exif");return this.interopOffset||(this.interopOffset=t.get(40965)),this.makerNote=t.get(37500),this.userComment=t.get(37510),this.options.sanitize&&(t.delete(40965),t.delete(37500),t.delete(37510)),this.unpack(t,41728),this.unpack(t,41729),t}unpack(t,e){let i=t.get(e);i&&i.length===1&&t.set(e,i[0])}async parseGpsBlock(){if(this.gps||(this.ifd0||await this.parseIfd0Block(),this.gpsOffset===void 0))return;let t=this.parseBlock(this.gpsOffset,"gps");return t&&t.has(2)&&t.has(4)&&(t.set("latitude",Ur(...t.get(2),t.get(1))),t.set("longitude",Ur(...t.get(4),t.get(3)))),t}async parseInteropBlock(){if(!this.interop&&(this.ifd0||await this.parseIfd0Block(),this.interopOffset!==void 0||this.exif||await this.parseExifBlock(),this.interopOffset!==void 0))return this.parseBlock(this.interopOffset,"interop")}async parseThumbnailBlock(t=!1){if(!this.ifd1&&!this.ifd1Parsed&&(!this.options.mergeOutput||t))return this.findIfd1Offset(),this.ifd1Offset>0&&(this.parseBlock(this.ifd1Offset,"ifd1"),this.ifd1Parsed=!0),this.ifd1}async extractThumbnail(){if(this.headerParsed||this.parseHeader(),this.ifd1Parsed||await this.parseThumbnailBlock(!0),this.ifd1===void 0)return;let t=this.ifd1.get(513),e=this.ifd1.get(514);return this.chunk.getUint8Array(t,e)}get image(){return this.ifd0}get thumbnail(){return this.ifd1}createOutput(){let t,e,i,r={};for(e of bt)if(t=this[e],!ms(t))if(i=this.canTranslate?this.translateBlock(t,e):Object.fromEntries(t),this.options.mergeOutput){if(e==="ifd1")continue;Object.assign(r,i)}else r[e]=i;return this.makerNote&&(r.makerNote=this.makerNote),this.userComment&&(r.userComment=this.userComment),r}assignToOutput(t,e){if(this.globalOptions.mergeOutput)Object.assign(t,e);else for(let[i,r]of Object.entries(e))this.assignObjectToOutput(t,i,r)}}function Ur(n,t,e,i){var r=n+t/60+e/3600;return i!=="S"&&i!=="W"||(r*=-1),r}J(Ei,"type","tiff"),J(Ei,"headerLength",10),Mt.set("tiff",Ei);var za=Object.freeze({__proto__:null,default:La,Exifr:Ne,fileParsers:ue,segmentParsers:Mt,fileReaders:de,tagKeys:Ct,tagValues:Jt,tagRevivers:Ie,createDictionary:wt,extendDictionary:an,fetchUrlAsArrayBuffer:rn,readBlobAsArrayBuffer:sn,chunkedProps:Pe,otherSegments:yn,segments:on,tiffBlocks:bt,segmentsAndBlocks:Ce,tiffExtractables:ke,inheritables:xn,allFormatters:Fe,Options:ln,parse:fi});const tr={ifd0:!1,ifd1:!1,exif:!1,gps:!1,interop:!1,sanitize:!1,reviveValues:!0,translateKeys:!1,translateValues:!1,mergeOutput:!1},er=Object.assign({},tr,{firstChunkSize:4e4,gps:[1,2,3,4]});async function xs(n){let t=new Ne(er);await t.read(n);let e=await t.parse();if(e&&e.gps){let{latitude:i,longitude:r}=e.gps;return{latitude:i,longitude:r}}}const nr=Object.assign({},tr,{tiff:!1,ifd1:!0,mergeOutput:!1});async function bs(n){let t=new Ne(nr);await t.read(n);let e=await t.extractThumbnail();return e&&di?hi.from(e):e}async function _s(n){let t=await this.thumbnail(n);if(t!==void 0){let e=new Blob([t]);return URL.createObjectURL(e)}}const ir=Object.assign({},tr,{firstChunkSize:4e4,ifd0:[274]});async function rr(n){let t=new Ne(ir);await t.read(n);let e=await t.parse();if(e&&e.ifd0)return e.ifd0[274]}const sr=Object.freeze({1:{dimensionSwapped:!1,scaleX:1,scaleY:1,deg:0,rad:0},2:{dimensionSwapped:!1,scaleX:-1,scaleY:1,deg:0,rad:0},3:{dimensionSwapped:!1,scaleX:1,scaleY:1,deg:180,rad:180*Math.PI/180},4:{dimensionSwapped:!1,scaleX:-1,scaleY:1,deg:180,rad:180*Math.PI/180},5:{dimensionSwapped:!0,scaleX:1,scaleY:-1,deg:90,rad:90*Math.PI/180},6:{dimensionSwapped:!0,scaleX:1,scaleY:1,deg:90,rad:90*Math.PI/180},7:{dimensionSwapped:!0,scaleX:1,scaleY:-1,deg:270,rad:270*Math.PI/180},8:{dimensionSwapped:!0,scaleX:1,scaleY:1,deg:270,rad:270*Math.PI/180}});let Je=!0,Ze=!0;if(typeof navigator=="object"){let n=navigator.userAgent;if(n.includes("iPad")||n.includes("iPhone")){let t=n.match(/OS (\d+)_(\d+)/);if(t){let[,e,i]=t;Je=Number(e)+.1*Number(i)<13.4,Ze=!1}}else if(n.includes("OS X 10")){let[,t]=n.match(/OS X 10[_.](\d+)/);Je=Ze=Number(t)<15}if(n.includes("Chrome/")){let[,t]=n.match(/Chrome\/(\d+)/);Je=Ze=Number(t)<81}else if(n.includes("Firefox/")){let[,t]=n.match(/Firefox\/(\d+)/);Je=Ze=Number(t)<77}}async function ws(n){let t=await rr(n);return Object.assign({canvas:Je,css:Ze},sr[t])}class Va extends Dt{constructor(...t){super(...t),J(this,"ranges",new Ga),this.byteLength!==0&&this.ranges.add(0,this.byteLength)}_tryExtend(t,e,i){if(t===0&&this.byteLength===0&&i){let r=new DataView(i.buffer||i,i.byteOffset,i.byteLength);this._swapDataView(r)}else{let r=t+e;if(r>this.byteLength){let{dataView:s}=this._extend(r);this._swapDataView(s)}}}_extend(t){let e;e=di?hi.allocUnsafe(t):new Uint8Array(t);let i=new DataView(e.buffer,e.byteOffset,e.byteLength);return e.set(new Uint8Array(this.buffer,this.byteOffset,this.byteLength),0),{uintView:e,dataView:i}}subarray(t,e,i=!1){return e=e||this._lengthToEnd(t),i&&this._tryExtend(t,e),this.ranges.add(t,e),super.subarray(t,e)}set(t,e,i=!1){i&&this._tryExtend(e,t.byteLength,t);let r=super.set(t,e);return this.ranges.add(e,r.byteLength),r}async ensureChunk(t,e){this.chunked&&(this.ranges.available(t,e)||await this.readChunk(t,e))}available(t,e){return this.ranges.available(t,e)}}class Ga{constructor(){J(this,"list",[])}get length(){return this.list.length}add(t,e,i=0){let r=t+e,s=this.list.filter(o=>Dr(t,o.offset,r)||Dr(t,o.end,r));if(s.length>0){t=Math.min(t,...s.map(a=>a.offset)),r=Math.max(r,...s.map(a=>a.end)),e=r-t;let o=s.shift();o.offset=t,o.length=e,o.end=r,this.list=this.list.filter(a=>!s.includes(a))}else this.list.push({offset:t,length:e,end:r})}available(t,e){let i=t+e;return this.list.some(r=>r.offset<=t&&i<=r.end)}}function Dr(n,t,e){return n<=t&&t<=e}class mi extends Va{constructor(t,e){super(0),J(this,"chunksRead",0),this.input=t,this.options=e}async readWhole(){this.chunked=!1,await this.readChunk(this.nextChunkOffset)}async readChunked(){this.chunked=!0,await this.readChunk(0,this.options.firstChunkSize)}async readNextChunk(t=this.nextChunkOffset){if(this.fullyRead)return this.chunksRead++,!1;let e=this.options.chunkSize,i=await this.readChunk(t,e);return!!i&&i.byteLength===e}async readChunk(t,e){if(this.chunksRead++,(e=this.safeWrapAddress(t,e))!==0)return this._readChunk(t,e)}safeWrapAddress(t,e){return this.size!==void 0&&t+e>this.size?Math.max(0,this.size-t):e}get nextChunkOffset(){if(this.ranges.list.length!==0)return this.ranges.list[0].length}get canReadNextChunk(){return this.chunksRead<this.options.chunkLimit}get fullyRead(){return this.size!==void 0&&this.nextChunkOffset===this.size}read(){return this.options.chunked?this.readChunked():this.readWhole()}close(){}}de.set("blob",class extends mi{async readWhole(){this.chunked=!1;let n=await sn(this.input);this._swapArrayBuffer(n)}readChunked(){return this.chunked=!0,this.size=this.input.size,super.readChunked()}async _readChunk(n,t){let e=t?n+t:void 0,i=this.input.slice(n,e),r=await sn(i);return this.set(r,n,!0)}});var Xa=Object.freeze({__proto__:null,default:za,Exifr:Ne,fileParsers:ue,segmentParsers:Mt,fileReaders:de,tagKeys:Ct,tagValues:Jt,tagRevivers:Ie,createDictionary:wt,extendDictionary:an,fetchUrlAsArrayBuffer:rn,readBlobAsArrayBuffer:sn,chunkedProps:Pe,otherSegments:yn,segments:on,tiffBlocks:bt,segmentsAndBlocks:Ce,tiffExtractables:ke,inheritables:xn,allFormatters:Fe,Options:ln,parse:fi,gpsOnlyOptions:er,gps:xs,thumbnailOnlyOptions:nr,thumbnail:bs,thumbnailUrl:_s,orientationOnlyOptions:ir,orientation:rr,rotations:sr,get rotateCanvas(){return Je},get rotateCss(){return Ze},rotation:ws});de.set("url",class extends mi{async readWhole(){this.chunked=!1;let n=await rn(this.input);n instanceof ArrayBuffer?this._swapArrayBuffer(n):n instanceof Uint8Array&&this._swapBuffer(n)}async _readChunk(n,t){let e=t?n+t-1:void 0,i=this.options.httpHeaders||{};(n||e)&&(i.range=`bytes=${[n,e].join("-")}`);let r=await Ji(this.input,{headers:i}),s=await r.arrayBuffer(),o=s.byteLength;if(r.status!==416)return o!==t&&(this.size=n+o),this.set(s,n,!0)}});Dt.prototype.getUint64=function(n){let t=this.getUint32(n),e=this.getUint32(n+4);return t<1048575?t<<32|e:typeof Wn!==void 0?(console.warn("Using BigInt because of type 64uint but JS can only handle 53b numbers."),Wn(t)<<Wn(32)|Wn(e)):void St("Trying to read 64b value but JS can only handle 53b numbers.")};class Ya extends pi{parseBoxes(t=0){let e=[];for(;t<this.file.byteLength-4;){let i=this.parseBoxHead(t);if(e.push(i),i.length===0)break;t+=i.length}return e}parseSubBoxes(t){t.boxes=this.parseBoxes(t.start)}findBox(t,e){return t.boxes===void 0&&this.parseSubBoxes(t),t.boxes.find(i=>i.kind===e)}parseBoxHead(t){let e=this.file.getUint32(t),i=this.file.getString(t+4,4),r=t+8;return e===1&&(e=this.file.getUint64(t+8),r+=8),{offset:t,length:e,kind:i,start:r}}parseBoxFullHead(t){if(t.version!==void 0)return;let e=this.file.getUint32(t.start);t.version=e>>24,t.start+=4}}class Ms extends Ya{static canHandle(t,e){if(e!==0)return!1;let i=t.getUint16(2);if(i>50)return!1;let r=16,s=[];for(;r<i;)s.push(t.getString(r,4)),r+=4;return s.includes(this.type)}async parse(){let t=this.file.getUint32(0),e=this.parseBoxHead(t);for(;e.kind!=="meta";)t+=e.length,await this.file.ensureChunk(t,16),e=this.parseBoxHead(t);await this.file.ensureChunk(e.offset,e.length),this.parseBoxFullHead(e),this.parseSubBoxes(e),this.options.icc.enabled&&await this.findIcc(e),this.options.tiff.enabled&&await this.findExif(e)}async registerSegment(t,e,i){await this.file.ensureChunk(e,i);let r=this.file.subarray(e,i);this.createParser(t,r)}async findIcc(t){let e=this.findBox(t,"iprp");if(e===void 0)return;let i=this.findBox(e,"ipco");if(i===void 0)return;let r=this.findBox(i,"colr");r!==void 0&&await this.registerSegment("icc",r.offset+12,r.length)}async findExif(t){let e=this.findBox(t,"iinf");if(e===void 0)return;let i=this.findBox(t,"iloc");if(i===void 0)return;let r=this.findExifLocIdInIinf(e),s=this.findExtentInIloc(i,r);if(s===void 0)return;let[o,a]=s;await this.file.ensureChunk(o,a);let l=4+this.file.getUint32(o);o+=l,a-=l,await this.registerSegment("tiff",o,a)}findExifLocIdInIinf(t){this.parseBoxFullHead(t);let e,i,r,s,o=t.start,a=this.file.getUint16(o);for(o+=2;a--;){if(e=this.parseBoxHead(o),this.parseBoxFullHead(e),i=e.start,e.version>=2&&(r=e.version===3?4:2,s=this.file.getString(i+r+2,4),s==="Exif"))return this.file.getUintBytes(i,r);o+=e.length}}get8bits(t){let e=this.file.getUint8(t);return[e>>4,15&e]}findExtentInIloc(t,e){this.parseBoxFullHead(t);let i=t.start,[r,s]=this.get8bits(i++),[o,a]=this.get8bits(i++),l=t.version===2?4:2,u=t.version===1||t.version===2?2:0,h=a+r+s,c=t.version===2?4:2,y=this.file.getUintBytes(i,c);for(i+=c;y--;){let f=this.file.getUintBytes(i,l);i+=l+u+2+o;let p=this.file.getUint16(i);if(i+=2,f===e)return p>1&&console.warn(`ILOC box has more than one extent but we're only processing one
Please create an issue at https://github.com/MikeKovarik/exifr with this file`),[this.file.getUintBytes(i+a,r),this.file.getUintBytes(i+a+r,s)];i+=p*h}}}class Ss extends Ms{}J(Ss,"type","heic");class Br extends Ms{}J(Br,"type","avif"),ue.set("heic",Ss),ue.set("avif",Br),wt(Ct,["ifd0","ifd1"],[[256,"ImageWidth"],[257,"ImageHeight"],[258,"BitsPerSample"],[259,"Compression"],[262,"PhotometricInterpretation"],[270,"ImageDescription"],[271,"Make"],[272,"Model"],[273,"StripOffsets"],[274,"Orientation"],[277,"SamplesPerPixel"],[278,"RowsPerStrip"],[279,"StripByteCounts"],[282,"XResolution"],[283,"YResolution"],[284,"PlanarConfiguration"],[296,"ResolutionUnit"],[301,"TransferFunction"],[305,"Software"],[306,"ModifyDate"],[315,"Artist"],[316,"HostComputer"],[317,"Predictor"],[318,"WhitePoint"],[319,"PrimaryChromaticities"],[513,"ThumbnailOffset"],[514,"ThumbnailLength"],[529,"YCbCrCoefficients"],[530,"YCbCrSubSampling"],[531,"YCbCrPositioning"],[532,"ReferenceBlackWhite"],[700,"ApplicationNotes"],[33432,"Copyright"],[33723,"IPTC"],[34665,"ExifIFD"],[34675,"ICC"],[34853,"GpsIFD"],[330,"SubIFD"],[40965,"InteropIFD"],[40091,"XPTitle"],[40092,"XPComment"],[40093,"XPAuthor"],[40094,"XPKeywords"],[40095,"XPSubject"]]),wt(Ct,"exif",[[33434,"ExposureTime"],[33437,"FNumber"],[34850,"ExposureProgram"],[34852,"SpectralSensitivity"],[34855,"ISO"],[34858,"TimeZoneOffset"],[34859,"SelfTimerMode"],[34864,"SensitivityType"],[34865,"StandardOutputSensitivity"],[34866,"RecommendedExposureIndex"],[34867,"ISOSpeed"],[34868,"ISOSpeedLatitudeyyy"],[34869,"ISOSpeedLatitudezzz"],[36864,"ExifVersion"],[36867,"DateTimeOriginal"],[36868,"CreateDate"],[36873,"GooglePlusUploadCode"],[36880,"OffsetTime"],[36881,"OffsetTimeOriginal"],[36882,"OffsetTimeDigitized"],[37121,"ComponentsConfiguration"],[37122,"CompressedBitsPerPixel"],[37377,"ShutterSpeedValue"],[37378,"ApertureValue"],[37379,"BrightnessValue"],[37380,"ExposureCompensation"],[37381,"MaxApertureValue"],[37382,"SubjectDistance"],[37383,"MeteringMode"],[37384,"LightSource"],[37385,"Flash"],[37386,"FocalLength"],[37393,"ImageNumber"],[37394,"SecurityClassification"],[37395,"ImageHistory"],[37396,"SubjectArea"],[37500,"MakerNote"],[37510,"UserComment"],[37520,"SubSecTime"],[37521,"SubSecTimeOriginal"],[37522,"SubSecTimeDigitized"],[37888,"AmbientTemperature"],[37889,"Humidity"],[37890,"Pressure"],[37891,"WaterDepth"],[37892,"Acceleration"],[37893,"CameraElevationAngle"],[40960,"FlashpixVersion"],[40961,"ColorSpace"],[40962,"ExifImageWidth"],[40963,"ExifImageHeight"],[40964,"RelatedSoundFile"],[41483,"FlashEnergy"],[41486,"FocalPlaneXResolution"],[41487,"FocalPlaneYResolution"],[41488,"FocalPlaneResolutionUnit"],[41492,"SubjectLocation"],[41493,"ExposureIndex"],[41495,"SensingMethod"],[41728,"FileSource"],[41729,"SceneType"],[41730,"CFAPattern"],[41985,"CustomRendered"],[41986,"ExposureMode"],[41987,"WhiteBalance"],[41988,"DigitalZoomRatio"],[41989,"FocalLengthIn35mmFormat"],[41990,"SceneCaptureType"],[41991,"GainControl"],[41992,"Contrast"],[41993,"Saturation"],[41994,"Sharpness"],[41996,"SubjectDistanceRange"],[42016,"ImageUniqueID"],[42032,"OwnerName"],[42033,"SerialNumber"],[42034,"LensInfo"],[42035,"LensMake"],[42036,"LensModel"],[42037,"LensSerialNumber"],[42080,"CompositeImage"],[42081,"CompositeImageCount"],[42082,"CompositeImageExposureTimes"],[42240,"Gamma"],[59932,"Padding"],[59933,"OffsetSchema"],[65e3,"OwnerName"],[65001,"SerialNumber"],[65002,"Lens"],[65100,"RawFile"],[65101,"Converter"],[65102,"WhiteBalance"],[65105,"Exposure"],[65106,"Shadows"],[65107,"Brightness"],[65108,"Contrast"],[65109,"Saturation"],[65110,"Sharpness"],[65111,"Smoothness"],[65112,"MoireFilter"],[40965,"InteropIFD"]]),wt(Ct,"gps",[[0,"GPSVersionID"],[1,"GPSLatitudeRef"],[2,"GPSLatitude"],[3,"GPSLongitudeRef"],[4,"GPSLongitude"],[5,"GPSAltitudeRef"],[6,"GPSAltitude"],[7,"GPSTimeStamp"],[8,"GPSSatellites"],[9,"GPSStatus"],[10,"GPSMeasureMode"],[11,"GPSDOP"],[12,"GPSSpeedRef"],[13,"GPSSpeed"],[14,"GPSTrackRef"],[15,"GPSTrack"],[16,"GPSImgDirectionRef"],[17,"GPSImgDirection"],[18,"GPSMapDatum"],[19,"GPSDestLatitudeRef"],[20,"GPSDestLatitude"],[21,"GPSDestLongitudeRef"],[22,"GPSDestLongitude"],[23,"GPSDestBearingRef"],[24,"GPSDestBearing"],[25,"GPSDestDistanceRef"],[26,"GPSDestDistance"],[27,"GPSProcessingMethod"],[28,"GPSAreaInformation"],[29,"GPSDateStamp"],[30,"GPSDifferential"],[31,"GPSHPositioningError"]]),wt(Jt,["ifd0","ifd1"],[[274,{1:"Horizontal (normal)",2:"Mirror horizontal",3:"Rotate 180",4:"Mirror vertical",5:"Mirror horizontal and rotate 270 CW",6:"Rotate 90 CW",7:"Mirror horizontal and rotate 90 CW",8:"Rotate 270 CW"}],[296,{1:"None",2:"inches",3:"cm"}]]);let In=wt(Jt,"exif",[[34850,{0:"Not defined",1:"Manual",2:"Normal program",3:"Aperture priority",4:"Shutter priority",5:"Creative program",6:"Action program",7:"Portrait mode",8:"Landscape mode"}],[37121,{0:"-",1:"Y",2:"Cb",3:"Cr",4:"R",5:"G",6:"B"}],[37383,{0:"Unknown",1:"Average",2:"CenterWeightedAverage",3:"Spot",4:"MultiSpot",5:"Pattern",6:"Partial",255:"Other"}],[37384,{0:"Unknown",1:"Daylight",2:"Fluorescent",3:"Tungsten (incandescent light)",4:"Flash",9:"Fine weather",10:"Cloudy weather",11:"Shade",12:"Daylight fluorescent (D 5700 - 7100K)",13:"Day white fluorescent (N 4600 - 5400K)",14:"Cool white fluorescent (W 3900 - 4500K)",15:"White fluorescent (WW 3200 - 3700K)",17:"Standard light A",18:"Standard light B",19:"Standard light C",20:"D55",21:"D65",22:"D75",23:"D50",24:"ISO studio tungsten",255:"Other"}],[37385,{0:"Flash did not fire",1:"Flash fired",5:"Strobe return light not detected",7:"Strobe return light detected",9:"Flash fired, compulsory flash mode",13:"Flash fired, compulsory flash mode, return light not detected",15:"Flash fired, compulsory flash mode, return light detected",16:"Flash did not fire, compulsory flash mode",24:"Flash did not fire, auto mode",25:"Flash fired, auto mode",29:"Flash fired, auto mode, return light not detected",31:"Flash fired, auto mode, return light detected",32:"No flash function",65:"Flash fired, red-eye reduction mode",69:"Flash fired, red-eye reduction mode, return light not detected",71:"Flash fired, red-eye reduction mode, return light detected",73:"Flash fired, compulsory flash mode, red-eye reduction mode",77:"Flash fired, compulsory flash mode, red-eye reduction mode, return light not detected",79:"Flash fired, compulsory flash mode, red-eye reduction mode, return light detected",89:"Flash fired, auto mode, red-eye reduction mode",93:"Flash fired, auto mode, return light not detected, red-eye reduction mode",95:"Flash fired, auto mode, return light detected, red-eye reduction mode"}],[41495,{1:"Not defined",2:"One-chip color area sensor",3:"Two-chip color area sensor",4:"Three-chip color area sensor",5:"Color sequential area sensor",7:"Trilinear sensor",8:"Color sequential linear sensor"}],[41728,{1:"Film Scanner",2:"Reflection Print Scanner",3:"Digital Camera"}],[41729,{1:"Directly photographed"}],[41985,{0:"Normal",1:"Custom",2:"HDR (no original saved)",3:"HDR (original saved)",4:"Original (for HDR)",6:"Panorama",7:"Portrait HDR",8:"Portrait"}],[41986,{0:"Auto",1:"Manual",2:"Auto bracket"}],[41987,{0:"Auto",1:"Manual"}],[41990,{0:"Standard",1:"Landscape",2:"Portrait",3:"Night",4:"Other"}],[41991,{0:"None",1:"Low gain up",2:"High gain up",3:"Low gain down",4:"High gain down"}],[41996,{0:"Unknown",1:"Macro",2:"Close",3:"Distant"}],[42080,{0:"Unknown",1:"Not a Composite Image",2:"General Composite Image",3:"Composite Image Captured While Shooting"}]]);const Or={1:"No absolute unit of measurement",2:"Inch",3:"Centimeter"};In.set(37392,Or),In.set(41488,Or);const Ui={0:"Normal",1:"Low",2:"High"};function zr(n){return typeof n=="object"&&n.length!==void 0?n[0]:n}function Vr(n){let t=Array.from(n).slice(1);return t[1]>15&&(t=t.map(e=>String.fromCharCode(e))),t[2]!=="0"&&t[2]!==0||t.pop(),t.join(".")}function Di(n){if(typeof n=="string"){var[t,e,i,r,s,o]=n.trim().split(/[-: ]/g).map(Number),a=new Date(t,e-1,i);return Number.isNaN(r)||Number.isNaN(s)||Number.isNaN(o)||(a.setHours(r),a.setMinutes(s),a.setSeconds(o)),Number.isNaN(+a)?n:a}}function Fn(n){if(typeof n=="string")return n;let t=[];if(n[1]===0&&n[n.length-1]===0)for(let e=0;e<n.length;e+=2)t.push(Gr(n[e+1],n[e]));else for(let e=0;e<n.length;e+=2)t.push(Gr(n[e],n[e+1]));return $e(String.fromCodePoint(...t))}function Gr(n,t){return n<<8|t}In.set(41992,Ui),In.set(41993,Ui),In.set(41994,Ui),wt(Ie,["ifd0","ifd1"],[[50827,function(n){return typeof n!="string"?gs(n):n}],[306,Di],[40091,Fn],[40092,Fn],[40093,Fn],[40094,Fn],[40095,Fn]]),wt(Ie,"exif",[[40960,Vr],[36864,Vr],[36867,Di],[36868,Di],[40962,zr],[40963,zr]]),wt(Ie,"gps",[[0,n=>Array.from(n).join(".")],[7,n=>Array.from(n).join(":")]]);class Bi extends $t{static canHandle(t,e){return t.getUint8(e+1)===225&&t.getUint32(e+4)===1752462448&&t.getString(e+4,20)==="http://ns.adobe.com/"}static headerLength(t,e){return t.getString(e+4,34)==="http://ns.adobe.com/xmp/extension/"?79:33}static findPosition(t,e){let i=super.findPosition(t,e);return i.multiSegment=i.extended=i.headerLength===79,i.multiSegment?(i.chunkCount=t.getUint8(e+72),i.chunkNumber=t.getUint8(e+76),t.getUint8(e+77)!==0&&i.chunkNumber++):(i.chunkCount=1/0,i.chunkNumber=-1),i}static handleMultiSegments(t){return t.map(e=>e.chunk.getString()).join("")}normalizeInput(t){return typeof t=="string"?t:Dt.from(t).getString()}parse(t=this.chunk){if(!this.localOptions.parse)return t;t=function(s){let o={},a={};for(let l of ks)o[l]=[],a[l]=0;return s.replace(qa,(l,u,h)=>{if(u==="<"){let c=++a[h];return o[h].push(c),`${l}#${c}`}return`${l}#${o[h].pop()}`})}(t);let e=pn.findAll(t,"rdf","Description");e.length===0&&e.push(new pn("rdf","Description",void 0,t));let i,r={};for(let s of e)for(let o of s.properties)i=Ha(o.ns,r),vs(o,i);return function(s){let o;for(let a in s)o=s[a]=ni(s[a]),o===void 0&&delete s[a];return ni(s)}(r)}assignToOutput(t,e){if(this.localOptions.parse)for(let[i,r]of Object.entries(e))switch(i){case"tiff":this.assignObjectToOutput(t,"ifd0",r);break;case"exif":this.assignObjectToOutput(t,"exif",r);break;case"xmlns":break;default:this.assignObjectToOutput(t,i,r)}else t.xmp=e}}J(Bi,"type","xmp"),J(Bi,"multiSegment",!0),Mt.set("xmp",Bi);class ri{static findAll(t){return Ps(t,/([a-zA-Z0-9-]+):([a-zA-Z0-9-]+)=("[^"]*"|'[^']*')/gm).map(ri.unpackMatch)}static unpackMatch(t){let e=t[1],i=t[2],r=t[3].slice(1,-1);return r=Cs(r),new ri(e,i,r)}constructor(t,e,i){this.ns=t,this.name=e,this.value=i}serialize(){return this.value}}class pn{static findAll(t,e,i){if(e!==void 0||i!==void 0){e=e||"[\\w\\d-]+",i=i||"[\\w\\d-]+";var r=new RegExp(`<(${e}):(${i})(#\\d+)?((\\s+?[\\w\\d-:]+=("[^"]*"|'[^']*'))*\\s*)(\\/>|>([\\s\\S]*?)<\\/\\1:\\2\\3>)`,"gm")}else r=/<([\w\d-]+):([\w\d-]+)(#\d+)?((\s+?[\w\d-:]+=("[^"]*"|'[^']*'))*\s*)(\/>|>([\s\S]*?)<\/\1:\2\3>)/gm;return Ps(t,r).map(pn.unpackMatch)}static unpackMatch(t){let e=t[1],i=t[2],r=t[4],s=t[8];return new pn(e,i,r,s)}constructor(t,e,i,r){this.ns=t,this.name=e,this.attrString=i,this.innerXml=r,this.attrs=ri.findAll(i),this.children=pn.findAll(r),this.value=this.children.length===0?Cs(r):void 0,this.properties=[...this.attrs,...this.children]}get isPrimitive(){return this.value!==void 0&&this.attrs.length===0&&this.children.length===0}get isListContainer(){return this.children.length===1&&this.children[0].isList}get isList(){let{ns:t,name:e}=this;return t==="rdf"&&(e==="Seq"||e==="Bag"||e==="Alt")}get isListItem(){return this.ns==="rdf"&&this.name==="li"}serialize(){if(this.properties.length===0&&this.value===void 0)return;if(this.isPrimitive)return this.value;if(this.isListContainer)return this.children[0].serialize();if(this.isList)return ja(this.children.map(Wa));if(this.isListItem&&this.children.length===1&&this.attrs.length===0)return this.children[0].serialize();let t={};for(let e of this.properties)vs(e,t);return this.value!==void 0&&(t.value=this.value),ni(t)}}function vs(n,t){let e=n.serialize();e!==void 0&&(t[n.name]=e)}var Wa=n=>n.serialize(),ja=n=>n.length===1?n[0]:n,Ha=(n,t)=>t[n]?t[n]:t[n]={};function Ps(n,t){let e,i=[];if(!n)return i;for(;(e=t.exec(n))!==null;)i.push(e);return i}function Cs(n){if(function(i){return i==null||i==="null"||i==="undefined"||i===""||i.trim()===""}(n))return;let t=Number(n);if(!Number.isNaN(t))return t;let e=n.toLowerCase();return e==="true"||e!=="false"&&n.trim()}const ks=["rdf:li","rdf:Seq","rdf:Bag","rdf:Alt","rdf:Description"],qa=new RegExp(`(<|\\/)(${ks.join("|")})`,"g");var Fs=Object.freeze({__proto__:null,default:Xa,Exifr:Ne,fileParsers:ue,segmentParsers:Mt,fileReaders:de,tagKeys:Ct,tagValues:Jt,tagRevivers:Ie,createDictionary:wt,extendDictionary:an,fetchUrlAsArrayBuffer:rn,readBlobAsArrayBuffer:sn,chunkedProps:Pe,otherSegments:yn,segments:on,tiffBlocks:bt,segmentsAndBlocks:Ce,tiffExtractables:ke,inheritables:xn,allFormatters:Fe,Options:ln,parse:fi,gpsOnlyOptions:er,gps:xs,thumbnailOnlyOptions:nr,thumbnail:bs,thumbnailUrl:_s,orientationOnlyOptions:ir,orientation:rr,rotations:sr,get rotateCanvas(){return Je},get rotateCss(){return Ze},rotation:ws});let Xr=ei("fs",n=>n.promises);de.set("fs",class extends mi{async readWhole(){this.chunked=!1,this.fs=await Xr;let n=await this.fs.readFile(this.input);this._swapBuffer(n)}async readChunked(){this.chunked=!0,this.fs=await Xr,await this.open(),await this.readChunk(0,this.options.firstChunkSize)}async open(){this.fh===void 0&&(this.fh=await this.fs.open(this.input,"r"),this.size=(await this.fh.stat(this.input)).size)}async _readChunk(n,t){this.fh===void 0&&await this.open(),n+t>this.size&&(t=this.size-n);var e=this.subarray(n,t,!0);return await this.fh.read(e.dataView,0,t,n),e}async close(){if(this.fh){let n=this.fh;this.fh=void 0,await n.close()}}});de.set("base64",class extends mi{constructor(...n){super(...n),this.input=this.input.replace(/^data:([^;]+);base64,/gim,""),this.size=this.input.length/4*3,this.input.endsWith("==")?this.size-=2:this.input.endsWith("=")&&(this.size-=1)}async _readChunk(n,t){let e,i,r=this.input;n===void 0?(n=0,e=0,i=0):(e=4*Math.floor(n/3),i=n-e/4*3),t===void 0&&(t=this.size);let s=n+t,o=e+4*Math.ceil(s/3);r=r.slice(e,o);let a=Math.min(t,this.size-n);if(di){let l=hi.from(r,"base64").slice(i,i+a);return this.set(l,n,!0)}{let l=this.subarray(n,a,!0),u=atob(r),h=l.toUint8();for(let c=0;c<a;c++)h[c]=u.charCodeAt(i+c);return l}}});class Yr extends pi{static canHandle(t,e){return e===18761||e===19789}extendOptions(t){let{ifd0:e,xmp:i,iptc:r,icc:s}=t;i.enabled&&e.deps.add(700),r.enabled&&e.deps.add(33723),s.enabled&&e.deps.add(34675),e.finalizeFilters()}async parse(){let{tiff:t,xmp:e,iptc:i,icc:r}=this.options;if(t.enabled||e.enabled||i.enabled||r.enabled){let s=Math.max(Xi(this.options),this.options.chunkSize);await this.file.ensureChunk(0,s),this.createParser("tiff",this.file),this.parsers.tiff.parseHeader(),await this.parsers.tiff.parseIfd0Block(),this.adaptTiffPropAsSegment("xmp"),this.adaptTiffPropAsSegment("iptc"),this.adaptTiffPropAsSegment("icc")}}adaptTiffPropAsSegment(t){if(this.parsers.tiff[t]){let e=this.parsers.tiff[t];this.injectSegment(t,e)}}}J(Yr,"type","tiff"),ue.set("tiff",Yr);let Qa=ei("zlib");const Ka=["ihdr","iccp","text","itxt","exif"];class Wr extends pi{constructor(...t){super(...t),J(this,"catchError",e=>this.errors.push(e)),J(this,"metaChunks",[]),J(this,"unknownChunks",[])}static canHandle(t,e){return e===35152&&t.getUint32(0)===2303741511&&t.getUint32(4)===218765834}async parse(){let{file:t}=this;await this.findPngChunksInRange(8,t.byteLength),await this.readSegments(this.metaChunks),this.findIhdr(),this.parseTextChunks(),await this.findExif().catch(this.catchError),await this.findXmp().catch(this.catchError),await this.findIcc().catch(this.catchError)}async findPngChunksInRange(t,e){let{file:i}=this;for(;t<e;){let r=i.getUint32(t),s=i.getUint32(t+4),o=i.getString(t+4,4).toLowerCase(),a=r+4+4+4,l={type:o,offset:t,length:a,start:t+4+4,size:r,marker:s};Ka.includes(o)?this.metaChunks.push(l):this.unknownChunks.push(l),t+=a}}parseTextChunks(){let t=this.metaChunks.filter(e=>e.type==="text");for(let e of t){let[i,r]=this.file.getString(e.start,e.size).split("\0");this.injectKeyValToIhdr(i,r)}}injectKeyValToIhdr(t,e){let i=this.parsers.ihdr;i&&i.raw.set(t,e)}findIhdr(){let t=this.metaChunks.find(e=>e.type==="ihdr");t&&this.options.ihdr.enabled!==!1&&this.createParser("ihdr",t.chunk)}async findExif(){let t=this.metaChunks.find(e=>e.type==="exif");t&&this.injectSegment("tiff",t.chunk)}async findXmp(){let t=this.metaChunks.filter(e=>e.type==="itxt");for(let e of t)e.chunk.getString(0,17)==="XML:com.adobe.xmp"&&this.injectSegment("xmp",e.chunk)}async findIcc(){let t=this.metaChunks.find(a=>a.type==="iccp");if(!t)return;let{chunk:e}=t,i=e.getUint8Array(0,81),r=0;for(;r<80&&i[r]!==0;)r++;let s=r+2,o=e.getString(0,r);if(this.injectKeyValToIhdr("ProfileName",o),ti){let a=await Qa,l=e.getUint8Array(s);l=a.inflateSync(l),this.injectSegment("icc",l)}}}J(Wr,"type","png"),ue.set("png",Wr),wt(Ct,"interop",[[1,"InteropIndex"],[2,"InteropVersion"],[4096,"RelatedImageFileFormat"],[4097,"RelatedImageWidth"],[4098,"RelatedImageHeight"]]),an(Ct,"ifd0",[[11,"ProcessingSoftware"],[254,"SubfileType"],[255,"OldSubfileType"],[263,"Thresholding"],[264,"CellWidth"],[265,"CellLength"],[266,"FillOrder"],[269,"DocumentName"],[280,"MinSampleValue"],[281,"MaxSampleValue"],[285,"PageName"],[286,"XPosition"],[287,"YPosition"],[290,"GrayResponseUnit"],[297,"PageNumber"],[321,"HalftoneHints"],[322,"TileWidth"],[323,"TileLength"],[332,"InkSet"],[337,"TargetPrinter"],[18246,"Rating"],[18249,"RatingPercent"],[33550,"PixelScale"],[34264,"ModelTransform"],[34377,"PhotoshopSettings"],[50706,"DNGVersion"],[50707,"DNGBackwardVersion"],[50708,"UniqueCameraModel"],[50709,"LocalizedCameraModel"],[50736,"DNGLensInfo"],[50739,"ShadowScale"],[50740,"DNGPrivateData"],[33920,"IntergraphMatrix"],[33922,"ModelTiePoint"],[34118,"SEMInfo"],[34735,"GeoTiffDirectory"],[34736,"GeoTiffDoubleParams"],[34737,"GeoTiffAsciiParams"],[50341,"PrintIM"],[50721,"ColorMatrix1"],[50722,"ColorMatrix2"],[50723,"CameraCalibration1"],[50724,"CameraCalibration2"],[50725,"ReductionMatrix1"],[50726,"ReductionMatrix2"],[50727,"AnalogBalance"],[50728,"AsShotNeutral"],[50729,"AsShotWhiteXY"],[50730,"BaselineExposure"],[50731,"BaselineNoise"],[50732,"BaselineSharpness"],[50734,"LinearResponseLimit"],[50735,"CameraSerialNumber"],[50741,"MakerNoteSafety"],[50778,"CalibrationIlluminant1"],[50779,"CalibrationIlluminant2"],[50781,"RawDataUniqueID"],[50827,"OriginalRawFileName"],[50828,"OriginalRawFileData"],[50831,"AsShotICCProfile"],[50832,"AsShotPreProfileMatrix"],[50833,"CurrentICCProfile"],[50834,"CurrentPreProfileMatrix"],[50879,"ColorimetricReference"],[50885,"SRawType"],[50898,"PanasonicTitle"],[50899,"PanasonicTitle2"],[50931,"CameraCalibrationSig"],[50932,"ProfileCalibrationSig"],[50933,"ProfileIFD"],[50934,"AsShotProfileName"],[50936,"ProfileName"],[50937,"ProfileHueSatMapDims"],[50938,"ProfileHueSatMapData1"],[50939,"ProfileHueSatMapData2"],[50940,"ProfileToneCurve"],[50941,"ProfileEmbedPolicy"],[50942,"ProfileCopyright"],[50964,"ForwardMatrix1"],[50965,"ForwardMatrix2"],[50966,"PreviewApplicationName"],[50967,"PreviewApplicationVersion"],[50968,"PreviewSettingsName"],[50969,"PreviewSettingsDigest"],[50970,"PreviewColorSpace"],[50971,"PreviewDateTime"],[50972,"RawImageDigest"],[50973,"OriginalRawFileDigest"],[50981,"ProfileLookTableDims"],[50982,"ProfileLookTableData"],[51043,"TimeCodes"],[51044,"FrameRate"],[51058,"TStop"],[51081,"ReelName"],[51089,"OriginalDefaultFinalSize"],[51090,"OriginalBestQualitySize"],[51091,"OriginalDefaultCropSize"],[51105,"CameraLabel"],[51107,"ProfileHueSatMapEncoding"],[51108,"ProfileLookTableEncoding"],[51109,"BaselineExposureOffset"],[51110,"DefaultBlackRender"],[51111,"NewRawImageDigest"],[51112,"RawToPreviewGain"]]);let jr=[[273,"StripOffsets"],[279,"StripByteCounts"],[288,"FreeOffsets"],[289,"FreeByteCounts"],[291,"GrayResponseCurve"],[292,"T4Options"],[293,"T6Options"],[300,"ColorResponseUnit"],[320,"ColorMap"],[324,"TileOffsets"],[325,"TileByteCounts"],[326,"BadFaxLines"],[327,"CleanFaxData"],[328,"ConsecutiveBadFaxLines"],[330,"SubIFD"],[333,"InkNames"],[334,"NumberofInks"],[336,"DotRange"],[338,"ExtraSamples"],[339,"SampleFormat"],[340,"SMinSampleValue"],[341,"SMaxSampleValue"],[342,"TransferRange"],[343,"ClipPath"],[344,"XClipPathUnits"],[345,"YClipPathUnits"],[346,"Indexed"],[347,"JPEGTables"],[351,"OPIProxy"],[400,"GlobalParametersIFD"],[401,"ProfileType"],[402,"FaxProfile"],[403,"CodingMethods"],[404,"VersionYear"],[405,"ModeNumber"],[433,"Decode"],[434,"DefaultImageColor"],[435,"T82Options"],[437,"JPEGTables"],[512,"JPEGProc"],[515,"JPEGRestartInterval"],[517,"JPEGLosslessPredictors"],[518,"JPEGPointTransforms"],[519,"JPEGQTables"],[520,"JPEGDCTables"],[521,"JPEGACTables"],[559,"StripRowCounts"],[999,"USPTOMiscellaneous"],[18247,"XP_DIP_XML"],[18248,"StitchInfo"],[28672,"SonyRawFileType"],[28688,"SonyToneCurve"],[28721,"VignettingCorrection"],[28722,"VignettingCorrParams"],[28724,"ChromaticAberrationCorrection"],[28725,"ChromaticAberrationCorrParams"],[28726,"DistortionCorrection"],[28727,"DistortionCorrParams"],[29895,"SonyCropTopLeft"],[29896,"SonyCropSize"],[32781,"ImageID"],[32931,"WangTag1"],[32932,"WangAnnotation"],[32933,"WangTag3"],[32934,"WangTag4"],[32953,"ImageReferencePoints"],[32954,"RegionXformTackPoint"],[32955,"WarpQuadrilateral"],[32956,"AffineTransformMat"],[32995,"Matteing"],[32996,"DataType"],[32997,"ImageDepth"],[32998,"TileDepth"],[33300,"ImageFullWidth"],[33301,"ImageFullHeight"],[33302,"TextureFormat"],[33303,"WrapModes"],[33304,"FovCot"],[33305,"MatrixWorldToScreen"],[33306,"MatrixWorldToCamera"],[33405,"Model2"],[33421,"CFARepeatPatternDim"],[33422,"CFAPattern2"],[33423,"BatteryLevel"],[33424,"KodakIFD"],[33445,"MDFileTag"],[33446,"MDScalePixel"],[33447,"MDColorTable"],[33448,"MDLabName"],[33449,"MDSampleInfo"],[33450,"MDPrepDate"],[33451,"MDPrepTime"],[33452,"MDFileUnits"],[33589,"AdventScale"],[33590,"AdventRevision"],[33628,"UIC1Tag"],[33629,"UIC2Tag"],[33630,"UIC3Tag"],[33631,"UIC4Tag"],[33918,"IntergraphPacketData"],[33919,"IntergraphFlagRegisters"],[33921,"INGRReserved"],[34016,"Site"],[34017,"ColorSequence"],[34018,"IT8Header"],[34019,"RasterPadding"],[34020,"BitsPerRunLength"],[34021,"BitsPerExtendedRunLength"],[34022,"ColorTable"],[34023,"ImageColorIndicator"],[34024,"BackgroundColorIndicator"],[34025,"ImageColorValue"],[34026,"BackgroundColorValue"],[34027,"PixelIntensityRange"],[34028,"TransparencyIndicator"],[34029,"ColorCharacterization"],[34030,"HCUsage"],[34031,"TrapIndicator"],[34032,"CMYKEquivalent"],[34152,"AFCP_IPTC"],[34232,"PixelMagicJBIGOptions"],[34263,"JPLCartoIFD"],[34306,"WB_GRGBLevels"],[34310,"LeafData"],[34687,"TIFF_FXExtensions"],[34688,"MultiProfiles"],[34689,"SharedData"],[34690,"T88Options"],[34732,"ImageLayer"],[34750,"JBIGOptions"],[34856,"Opto-ElectricConvFactor"],[34857,"Interlace"],[34908,"FaxRecvParams"],[34909,"FaxSubAddress"],[34910,"FaxRecvTime"],[34929,"FedexEDR"],[34954,"LeafSubIFD"],[37387,"FlashEnergy"],[37388,"SpatialFrequencyResponse"],[37389,"Noise"],[37390,"FocalPlaneXResolution"],[37391,"FocalPlaneYResolution"],[37392,"FocalPlaneResolutionUnit"],[37397,"ExposureIndex"],[37398,"TIFF-EPStandardID"],[37399,"SensingMethod"],[37434,"CIP3DataFile"],[37435,"CIP3Sheet"],[37436,"CIP3Side"],[37439,"StoNits"],[37679,"MSDocumentText"],[37680,"MSPropertySetStorage"],[37681,"MSDocumentTextPosition"],[37724,"ImageSourceData"],[40965,"InteropIFD"],[40976,"SamsungRawPointersOffset"],[40977,"SamsungRawPointersLength"],[41217,"SamsungRawByteOrder"],[41218,"SamsungRawUnknown"],[41484,"SpatialFrequencyResponse"],[41485,"Noise"],[41489,"ImageNumber"],[41490,"SecurityClassification"],[41491,"ImageHistory"],[41494,"TIFF-EPStandardID"],[41995,"DeviceSettingDescription"],[42112,"GDALMetadata"],[42113,"GDALNoData"],[44992,"ExpandSoftware"],[44993,"ExpandLens"],[44994,"ExpandFilm"],[44995,"ExpandFilterLens"],[44996,"ExpandScanner"],[44997,"ExpandFlashLamp"],[46275,"HasselbladRawImage"],[48129,"PixelFormat"],[48130,"Transformation"],[48131,"Uncompressed"],[48132,"ImageType"],[48256,"ImageWidth"],[48257,"ImageHeight"],[48258,"WidthResolution"],[48259,"HeightResolution"],[48320,"ImageOffset"],[48321,"ImageByteCount"],[48322,"AlphaOffset"],[48323,"AlphaByteCount"],[48324,"ImageDataDiscard"],[48325,"AlphaDataDiscard"],[50215,"OceScanjobDesc"],[50216,"OceApplicationSelector"],[50217,"OceIDNumber"],[50218,"OceImageLogic"],[50255,"Annotations"],[50459,"HasselbladExif"],[50547,"OriginalFileName"],[50560,"USPTOOriginalContentType"],[50656,"CR2CFAPattern"],[50710,"CFAPlaneColor"],[50711,"CFALayout"],[50712,"LinearizationTable"],[50713,"BlackLevelRepeatDim"],[50714,"BlackLevel"],[50715,"BlackLevelDeltaH"],[50716,"BlackLevelDeltaV"],[50717,"WhiteLevel"],[50718,"DefaultScale"],[50719,"DefaultCropOrigin"],[50720,"DefaultCropSize"],[50733,"BayerGreenSplit"],[50737,"ChromaBlurRadius"],[50738,"AntiAliasStrength"],[50752,"RawImageSegmentation"],[50780,"BestQualityScale"],[50784,"AliasLayerMetadata"],[50829,"ActiveArea"],[50830,"MaskedAreas"],[50935,"NoiseReductionApplied"],[50974,"SubTileBlockSize"],[50975,"RowInterleaveFactor"],[51008,"OpcodeList1"],[51009,"OpcodeList2"],[51022,"OpcodeList3"],[51041,"NoiseProfile"],[51114,"CacheVersion"],[51125,"DefaultUserCrop"],[51157,"NikonNEFInfo"],[65024,"KdcIFD"]];an(Ct,"ifd0",jr),an(Ct,"exif",jr),wt(Jt,"gps",[[23,{M:"Magnetic North",T:"True North"}],[25,{K:"Kilometers",M:"Miles",N:"Nautical Miles"}]]);class Oi extends $t{static canHandle(t,e){return t.getUint8(e+1)===224&&t.getUint32(e+4)===1246120262&&t.getUint8(e+8)===0}parse(){return this.parseTags(),this.translate(),this.output}parseTags(){this.raw=new Map([[0,this.chunk.getUint16(0)],[2,this.chunk.getUint8(2)],[3,this.chunk.getUint16(3)],[5,this.chunk.getUint16(5)],[7,this.chunk.getUint8(7)],[8,this.chunk.getUint8(8)]])}}J(Oi,"type","jfif"),J(Oi,"headerLength",9),Mt.set("jfif",Oi),wt(Ct,"jfif",[[0,"JFIFVersion"],[2,"ResolutionUnit"],[3,"XResolution"],[5,"YResolution"],[7,"ThumbnailWidth"],[8,"ThumbnailHeight"]]);class Hr extends $t{parse(){return this.parseTags(),this.translate(),this.output}parseTags(){this.raw=new Map([[0,this.chunk.getUint32(0)],[4,this.chunk.getUint32(4)],[8,this.chunk.getUint8(8)],[9,this.chunk.getUint8(9)],[10,this.chunk.getUint8(10)],[11,this.chunk.getUint8(11)],[12,this.chunk.getUint8(12)],...Array.from(this.raw)])}}J(Hr,"type","ihdr"),Mt.set("ihdr",Hr),wt(Ct,"ihdr",[[0,"ImageWidth"],[4,"ImageHeight"],[8,"BitDepth"],[9,"ColorType"],[10,"Compression"],[11,"Filter"],[12,"Interlace"]]),wt(Jt,"ihdr",[[9,{0:"Grayscale",2:"RGB",3:"Palette",4:"Grayscale with Alpha",6:"RGB with Alpha",DEFAULT:"Unknown"}],[10,{0:"Deflate/Inflate",DEFAULT:"Unknown"}],[11,{0:"Adaptive",DEFAULT:"Unknown"}],[12,{0:"Noninterlaced",1:"Adam7 Interlace",DEFAULT:"Unknown"}]]);class Zn extends $t{static canHandle(t,e){return t.getUint8(e+1)===226&&t.getUint32(e+4)===1229144927}static findPosition(t,e){let i=super.findPosition(t,e);return i.chunkNumber=t.getUint8(e+16),i.chunkCount=t.getUint8(e+17),i.multiSegment=i.chunkCount>1,i}static handleMultiSegments(t){return function(e){let i=function(r){let s=r[0].constructor,o=0;for(let u of r)o+=u.length;let a=new s(o),l=0;for(let u of r)a.set(u,l),l+=u.length;return a}(e.map(r=>r.chunk.toUint8()));return new Dt(i)}(t)}parse(){return this.raw=new Map,this.parseHeader(),this.parseTags(),this.translate(),this.output}parseHeader(){let{raw:t}=this;this.chunk.byteLength<84&&St("ICC header is too short");for(let[e,i]of Object.entries($a)){e=parseInt(e,10);let r=i(this.chunk,e);r!=="\0\0\0\0"&&t.set(e,r)}}parseTags(){let t,e,i,r,s,{raw:o}=this,a=this.chunk.getUint32(128),l=132,u=this.chunk.byteLength;for(;a--;){if(t=this.chunk.getString(l,4),e=this.chunk.getUint32(l+4),i=this.chunk.getUint32(l+8),r=this.chunk.getString(e,4),e+i>u)return void console.warn("reached the end of the first ICC chunk. Enable options.tiff.multiSegment to read all ICC segments.");s=this.parseTag(r,e,i),s!==void 0&&s!=="\0\0\0\0"&&o.set(t,s),l+=12}}parseTag(t,e,i){switch(t){case"desc":return this.parseDesc(e);case"mluc":return this.parseMluc(e);case"text":return this.parseText(e,i);case"sig ":return this.parseSig(e)}if(!(e+i>this.chunk.byteLength))return this.chunk.getUint8Array(e,i)}parseDesc(t){let e=this.chunk.getUint32(t+8)-1;return $e(this.chunk.getString(t+12,e))}parseText(t,e){return $e(this.chunk.getString(t+8,e-8))}parseSig(t){return $e(this.chunk.getString(t+8,4))}parseMluc(t){let{chunk:e}=this,i=e.getUint32(t+8),r=e.getUint32(t+12),s=t+16,o=[];for(let a=0;a<i;a++){let l=e.getString(s+0,2),u=e.getString(s+2,2),h=e.getUint32(s+4),c=e.getUint32(s+8)+t,y=$e(e.getUnicodeString(c,h));o.push({lang:l,country:u,text:y}),s+=r}return i===1?o[0].text:o}translateValue(t,e){return typeof t=="string"?e[t]||e[t.toLowerCase()]||t:e[t]||t}}J(Zn,"type","icc"),J(Zn,"multiSegment",!0),J(Zn,"headerLength",18);const $a={4:ye,8:function(n,t){return[n.getUint8(t),n.getUint8(t+1)>>4,n.getUint8(t+1)%16].map(e=>e.toString(10)).join(".")},12:ye,16:ye,20:ye,24:function(n,t){const e=n.getUint16(t),i=n.getUint16(t+2)-1,r=n.getUint16(t+4),s=n.getUint16(t+6),o=n.getUint16(t+8),a=n.getUint16(t+10);return new Date(Date.UTC(e,i,r,s,o,a))},36:ye,40:ye,48:ye,52:ye,64:(n,t)=>n.getUint32(t),80:ye};function ye(n,t){return $e(n.getString(t,4))}Mt.set("icc",Zn),wt(Ct,"icc",[[4,"ProfileCMMType"],[8,"ProfileVersion"],[12,"ProfileClass"],[16,"ColorSpaceData"],[20,"ProfileConnectionSpace"],[24,"ProfileDateTime"],[36,"ProfileFileSignature"],[40,"PrimaryPlatform"],[44,"CMMFlags"],[48,"DeviceManufacturer"],[52,"DeviceModel"],[56,"DeviceAttributes"],[64,"RenderingIntent"],[68,"ConnectionSpaceIlluminant"],[80,"ProfileCreator"],[84,"ProfileID"],["Header","ProfileHeader"],["MS00","WCSProfiles"],["bTRC","BlueTRC"],["bXYZ","BlueMatrixColumn"],["bfd","UCRBG"],["bkpt","MediaBlackPoint"],["calt","CalibrationDateTime"],["chad","ChromaticAdaptation"],["chrm","Chromaticity"],["ciis","ColorimetricIntentImageState"],["clot","ColorantTableOut"],["clro","ColorantOrder"],["clrt","ColorantTable"],["cprt","ProfileCopyright"],["crdi","CRDInfo"],["desc","ProfileDescription"],["devs","DeviceSettings"],["dmdd","DeviceModelDesc"],["dmnd","DeviceMfgDesc"],["dscm","ProfileDescriptionML"],["fpce","FocalPlaneColorimetryEstimates"],["gTRC","GreenTRC"],["gXYZ","GreenMatrixColumn"],["gamt","Gamut"],["kTRC","GrayTRC"],["lumi","Luminance"],["meas","Measurement"],["meta","Metadata"],["mmod","MakeAndModel"],["ncl2","NamedColor2"],["ncol","NamedColor"],["ndin","NativeDisplayInfo"],["pre0","Preview0"],["pre1","Preview1"],["pre2","Preview2"],["ps2i","PS2RenderingIntent"],["ps2s","PostScript2CSA"],["psd0","PostScript2CRD0"],["psd1","PostScript2CRD1"],["psd2","PostScript2CRD2"],["psd3","PostScript2CRD3"],["pseq","ProfileSequenceDesc"],["psid","ProfileSequenceIdentifier"],["psvm","PS2CRDVMSize"],["rTRC","RedTRC"],["rXYZ","RedMatrixColumn"],["resp","OutputResponse"],["rhoc","ReflectionHardcopyOrigColorimetry"],["rig0","PerceptualRenderingIntentGamut"],["rig2","SaturationRenderingIntentGamut"],["rpoc","ReflectionPrintOutputColorimetry"],["sape","SceneAppearanceEstimates"],["scoe","SceneColorimetryEstimates"],["scrd","ScreeningDesc"],["scrn","Screening"],["targ","CharTarget"],["tech","Technology"],["vcgt","VideoCardGamma"],["view","ViewingConditions"],["vued","ViewingCondDesc"],["wtpt","MediaWhitePoint"]]);const jn={"4d2p":"Erdt Systems",AAMA:"Aamazing Technologies",ACER:"Acer",ACLT:"Acolyte Color Research",ACTI:"Actix Sytems",ADAR:"Adara Technology",ADBE:"Adobe",ADI:"ADI Systems",AGFA:"Agfa Graphics",ALMD:"Alps Electric",ALPS:"Alps Electric",ALWN:"Alwan Color Expertise",AMTI:"Amiable Technologies",AOC:"AOC International",APAG:"Apago",APPL:"Apple Computer",AST:"AST","AT&T":"AT&T",BAEL:"BARBIERI electronic",BRCO:"Barco NV",BRKP:"Breakpoint",BROT:"Brother",BULL:"Bull",BUS:"Bus Computer Systems","C-IT":"C-Itoh",CAMR:"Intel",CANO:"Canon",CARR:"Carroll Touch",CASI:"Casio",CBUS:"Colorbus PL",CEL:"Crossfield",CELx:"Crossfield",CGS:"CGS Publishing Technologies International",CHM:"Rochester Robotics",CIGL:"Colour Imaging Group, London",CITI:"Citizen",CL00:"Candela",CLIQ:"Color IQ",CMCO:"Chromaco",CMiX:"CHROMiX",COLO:"Colorgraphic Communications",COMP:"Compaq",COMp:"Compeq/Focus Technology",CONR:"Conrac Display Products",CORD:"Cordata Technologies",CPQ:"Compaq",CPRO:"ColorPro",CRN:"Cornerstone",CTX:"CTX International",CVIS:"ColorVision",CWC:"Fujitsu Laboratories",DARI:"Darius Technology",DATA:"Dataproducts",DCP:"Dry Creek Photo",DCRC:"Digital Contents Resource Center, Chung-Ang University",DELL:"Dell Computer",DIC:"Dainippon Ink and Chemicals",DICO:"Diconix",DIGI:"Digital","DL&C":"Digital Light & Color",DPLG:"Doppelganger",DS:"Dainippon Screen",DSOL:"DOOSOL",DUPN:"DuPont",EPSO:"Epson",ESKO:"Esko-Graphics",ETRI:"Electronics and Telecommunications Research Institute",EVER:"Everex Systems",EXAC:"ExactCODE",Eizo:"Eizo",FALC:"Falco Data Products",FF:"Fuji Photo Film",FFEI:"FujiFilm Electronic Imaging",FNRD:"Fnord Software",FORA:"Fora",FORE:"Forefront Technology",FP:"Fujitsu",FPA:"WayTech Development",FUJI:"Fujitsu",FX:"Fuji Xerox",GCC:"GCC Technologies",GGSL:"Global Graphics Software",GMB:"Gretagmacbeth",GMG:"GMG",GOLD:"GoldStar Technology",GOOG:"Google",GPRT:"Giantprint",GTMB:"Gretagmacbeth",GVC:"WayTech Development",GW2K:"Sony",HCI:"HCI",HDM:"Heidelberger Druckmaschinen",HERM:"Hermes",HITA:"Hitachi America",HP:"Hewlett-Packard",HTC:"Hitachi",HiTi:"HiTi Digital",IBM:"IBM",IDNT:"Scitex",IEC:"Hewlett-Packard",IIYA:"Iiyama North America",IKEG:"Ikegami Electronics",IMAG:"Image Systems",IMI:"Ingram Micro",INTC:"Intel",INTL:"N/A (INTL)",INTR:"Intra Electronics",IOCO:"Iocomm International Technology",IPS:"InfoPrint Solutions Company",IRIS:"Scitex",ISL:"Ichikawa Soft Laboratory",ITNL:"N/A (ITNL)",IVM:"IVM",IWAT:"Iwatsu Electric",Idnt:"Scitex",Inca:"Inca Digital Printers",Iris:"Scitex",JPEG:"Joint Photographic Experts Group",JSFT:"Jetsoft Development",JVC:"JVC Information Products",KART:"Scitex",KFC:"KFC Computek Components",KLH:"KLH Computers",KMHD:"Konica Minolta",KNCA:"Konica",KODA:"Kodak",KYOC:"Kyocera",Kart:"Scitex",LCAG:"Leica",LCCD:"Leeds Colour",LDAK:"Left Dakota",LEAD:"Leading Technology",LEXM:"Lexmark International",LINK:"Link Computer",LINO:"Linotronic",LITE:"Lite-On",Leaf:"Leaf",Lino:"Linotronic",MAGC:"Mag Computronic",MAGI:"MAG Innovision",MANN:"Mannesmann",MICN:"Micron Technology",MICR:"Microtek",MICV:"Microvitec",MINO:"Minolta",MITS:"Mitsubishi Electronics America",MITs:"Mitsuba",MNLT:"Minolta",MODG:"Modgraph",MONI:"Monitronix",MONS:"Monaco Systems",MORS:"Morse Technology",MOTI:"Motive Systems",MSFT:"Microsoft",MUTO:"MUTOH INDUSTRIES",Mits:"Mitsubishi Electric",NANA:"NANAO",NEC:"NEC",NEXP:"NexPress Solutions",NISS:"Nissei Sangyo America",NKON:"Nikon",NONE:"none",OCE:"Oce Technologies",OCEC:"OceColor",OKI:"Oki",OKID:"Okidata",OKIP:"Okidata",OLIV:"Olivetti",OLYM:"Olympus",ONYX:"Onyx Graphics",OPTI:"Optiquest",PACK:"Packard Bell",PANA:"Matsushita Electric Industrial",PANT:"Pantone",PBN:"Packard Bell",PFU:"PFU",PHIL:"Philips Consumer Electronics",PNTX:"HOYA",POne:"Phase One A/S",PREM:"Premier Computer Innovations",PRIN:"Princeton Graphic Systems",PRIP:"Princeton Publishing Labs",QLUX:"Hong Kong",QMS:"QMS",QPCD:"QPcard AB",QUAD:"QuadLaser",QUME:"Qume",RADI:"Radius",RDDx:"Integrated Color Solutions",RDG:"Roland DG",REDM:"REDMS Group",RELI:"Relisys",RGMS:"Rolf Gierling Multitools",RICO:"Ricoh",RNLD:"Edmund Ronald",ROYA:"Royal",RPC:"Ricoh Printing Systems",RTL:"Royal Information Electronics",SAMP:"Sampo",SAMS:"Samsung",SANT:"Jaime Santana Pomares",SCIT:"Scitex",SCRN:"Dainippon Screen",SDP:"Scitex",SEC:"Samsung",SEIK:"Seiko Instruments",SEIk:"Seikosha",SGUY:"ScanGuy.com",SHAR:"Sharp Laboratories",SICC:"International Color Consortium",SONY:"Sony",SPCL:"SpectraCal",STAR:"Star",STC:"Sampo Technology",Scit:"Scitex",Sdp:"Scitex",Sony:"Sony",TALO:"Talon Technology",TAND:"Tandy",TATU:"Tatung",TAXA:"TAXAN America",TDS:"Tokyo Denshi Sekei",TECO:"TECO Information Systems",TEGR:"Tegra",TEKT:"Tektronix",TI:"Texas Instruments",TMKR:"TypeMaker",TOSB:"Toshiba",TOSH:"Toshiba",TOTK:"TOTOKU ELECTRIC",TRIU:"Triumph",TSBT:"Toshiba",TTX:"TTX Computer Products",TVM:"TVM Professional Monitor",TW:"TW Casper",ULSX:"Ulead Systems",UNIS:"Unisys",UTZF:"Utz Fehlau & Sohn",VARI:"Varityper",VIEW:"Viewsonic",VISL:"Visual communication",VIVO:"Vivo Mobile Communication",WANG:"Wang",WLBR:"Wilbur Imaging",WTG2:"Ware To Go",WYSE:"WYSE Technology",XERX:"Xerox",XRIT:"X-Rite",ZRAN:"Zoran",Zebr:"Zebra Technologies",appl:"Apple Computer",bICC:"basICColor",berg:"bergdesign",ceyd:"Integrated Color Solutions",clsp:"MacDermid ColorSpan",ds:"Dainippon Screen",dupn:"DuPont",ffei:"FujiFilm Electronic Imaging",flux:"FluxData",iris:"Scitex",kart:"Scitex",lcms:"Little CMS",lino:"Linotronic",none:"none",ob4d:"Erdt Systems",obic:"Medigraph",quby:"Qubyx Sarl",scit:"Scitex",scrn:"Dainippon Screen",sdp:"Scitex",siwi:"SIWI GRAFIKA",yxym:"YxyMaster"},qr={scnr:"Scanner",mntr:"Monitor",prtr:"Printer",link:"Device Link",abst:"Abstract",spac:"Color Space Conversion Profile",nmcl:"Named Color",cenc:"ColorEncodingSpace profile",mid:"MultiplexIdentification profile",mlnk:"MultiplexLink profile",mvis:"MultiplexVisualization profile",nkpf:"Nikon Input Device Profile (NON-STANDARD!)"};wt(Jt,"icc",[[4,jn],[12,qr],[40,Object.assign({},jn,qr)],[48,jn],[80,jn],[64,{0:"Perceptual",1:"Relative Colorimetric",2:"Saturation",3:"Absolute Colorimetric"}],["tech",{amd:"Active Matrix Display",crt:"Cathode Ray Tube Display",kpcd:"Photo CD",pmd:"Passive Matrix Display",dcam:"Digital Camera",dcpj:"Digital Cinema Projector",dmpc:"Digital Motion Picture Camera",dsub:"Dye Sublimation Printer",epho:"Electrophotographic Printer",esta:"Electrostatic Printer",flex:"Flexography",fprn:"Film Writer",fscn:"Film Scanner",grav:"Gravure",ijet:"Ink Jet Printer",imgs:"Photo Image Setter",mpfr:"Motion Picture Film Recorder",mpfs:"Motion Picture Film Scanner",offs:"Offset Lithography",pjtv:"Projection Television",rpho:"Photographic Paper Printer",rscn:"Reflective Scanner",silk:"Silkscreen",twax:"Thermal Wax Printer",vidc:"Video Camera",vidm:"Video Monitor"}]]);class Hn extends $t{static canHandle(t,e,i){return t.getUint8(e+1)===237&&t.getString(e+4,9)==="Photoshop"&&this.containsIptc8bim(t,e,i)!==void 0}static headerLength(t,e,i){let r,s=this.containsIptc8bim(t,e,i);if(s!==void 0)return r=t.getUint8(e+s+7),r%2!=0&&(r+=1),r===0&&(r=4),s+8+r}static containsIptc8bim(t,e,i){for(let r=0;r<i;r++)if(this.isIptcSegmentHead(t,e+r))return r}static isIptcSegmentHead(t,e){return t.getUint8(e)===56&&t.getUint32(e)===943868237&&t.getUint16(e+4)===1028}parse(){let{raw:t}=this,e=this.chunk.byteLength-1,i=!1;for(let r=0;r<e;r++)if(this.chunk.getUint8(r)===28&&this.chunk.getUint8(r+1)===2){i=!0;let s=this.chunk.getUint16(r+3),o=this.chunk.getUint8(r+2),a=this.chunk.getLatin1String(r+5,s);t.set(o,this.pluralizeValue(t.get(o),a)),r+=4+s}else if(i)break;return this.translate(),this.output}pluralizeValue(t,e){return t!==void 0?t instanceof Array?(t.push(e),t):[t,e]:e}}J(Hn,"type","iptc"),J(Hn,"translateValues",!1),J(Hn,"reviveValues",!1),Mt.set("iptc",Hn),wt(Ct,"iptc",[[0,"ApplicationRecordVersion"],[3,"ObjectTypeReference"],[4,"ObjectAttributeReference"],[5,"ObjectName"],[7,"EditStatus"],[8,"EditorialUpdate"],[10,"Urgency"],[12,"SubjectReference"],[15,"Category"],[20,"SupplementalCategories"],[22,"FixtureIdentifier"],[25,"Keywords"],[26,"ContentLocationCode"],[27,"ContentLocationName"],[30,"ReleaseDate"],[35,"ReleaseTime"],[37,"ExpirationDate"],[38,"ExpirationTime"],[40,"SpecialInstructions"],[42,"ActionAdvised"],[45,"ReferenceService"],[47,"ReferenceDate"],[50,"ReferenceNumber"],[55,"DateCreated"],[60,"TimeCreated"],[62,"DigitalCreationDate"],[63,"DigitalCreationTime"],[65,"OriginatingProgram"],[70,"ProgramVersion"],[75,"ObjectCycle"],[80,"Byline"],[85,"BylineTitle"],[90,"City"],[92,"Sublocation"],[95,"State"],[100,"CountryCode"],[101,"Country"],[103,"OriginalTransmissionReference"],[105,"Headline"],[110,"Credit"],[115,"Source"],[116,"CopyrightNotice"],[118,"Contact"],[120,"Caption"],[121,"LocalCaption"],[122,"Writer"],[125,"RasterizedCaption"],[130,"ImageType"],[131,"ImageOrientation"],[135,"LanguageIdentifier"],[150,"AudioType"],[151,"AudioSamplingRate"],[152,"AudioSamplingResolution"],[153,"AudioDuration"],[154,"AudioOutcue"],[184,"JobID"],[185,"MasterDocumentID"],[186,"ShortDocumentID"],[187,"UniqueDocumentID"],[188,"OwnerID"],[200,"ObjectPreviewFileFormat"],[201,"ObjectPreviewFileVersion"],[202,"ObjectPreviewData"],[221,"Prefs"],[225,"ClassifyState"],[228,"SimilarityIndex"],[230,"DocumentNotes"],[231,"DocumentHistory"],[232,"ExifCameraInfo"],[255,"CatalogSets"]]),wt(Jt,"iptc",[[10,{0:"0 (reserved)",1:"1 (most urgent)",2:"2",3:"3",4:"4",5:"5 (normal urgency)",6:"6",7:"7",8:"8 (least urgent)",9:"9 (user-defined priority)"}],[75,{a:"Morning",b:"Both Morning and Evening",p:"Evening"}],[131,{L:"Landscape",P:"Portrait",S:"Square"}]]);function Ja(n,t,e){const i=e.type===1||e.type===2||e.type===7?1:e.type===3?2:e.type===4||e.type===9||e.type===11||e.type===13?4:e.type===5||e.type===10||e.type===12?8:0;if(i===0)throw new Error(`Unsupported TIFF field type: ${e.type}`);const s=e.count*i<=4?e.valueFieldOffset:n.getUint32(e.valueFieldOffset,t),o=[];for(let a=0;a<e.count;a++){const l=s+a*i;if(l<0||l+i>n.byteLength)throw new Error("Invalid TIFF field offset");if(e.type===1)o.push(n.getUint8(l));else if(e.type===2||e.type===7)o.push(n.getUint8(l));else if(e.type===3)o.push(n.getUint16(l,t));else if(e.type===5){const u=n.getUint32(l+4,t);o.push(u?n.getUint32(l,t)/u:0)}else if(e.type===13)o.push(n.getUint32(l,t));else if(e.type===9)o.push(n.getInt32(l,t));else if(e.type===10){const u=n.getInt32(l+4,t);o.push(u?n.getInt32(l,t)/u:0)}else e.type===11?o.push(n.getFloat32(l,t)):e.type===12?o.push(n.getFloat64(l,t)):o.push(n.getUint32(l,t))}return o}function Qr(n,t,e){if(!Number.isFinite(e)||e<=0||e+2>n.byteLength)throw new Error("Invalid TIFF IFD offset");const i=n.getUint16(e,t);if(e+2+i*12+4>n.byteLength)throw new Error("Corrupt TIFF IFD");const s=new Map;for(let o=0;o<i;o++){const a=e+2+o*12,l=n.getUint16(a,t),u=n.getUint16(a+2,t),h=n.getUint32(a+4,t);s.set(l,{type:u,count:h,valueFieldOffset:a+8})}return s}function Ut(n,t,e,i){const r=e.get(i);return r?Ja(n,t,r):[]}const As=(...n)=>{for(const t of n){if(Array.isArray(t)||ArrayBuffer.isView(t)){const i=Array.from(t).map(Number).filter(r=>Number.isFinite(r)&&r>0);if(i.length)return Math.max(...i);continue}const e=Number(t);if(Number.isFinite(e)&&e>0)return e}return null},Za=n=>{const e=(Array.isArray(n)||ArrayBuffer.isView(n)?Array.from(n).map(Number):[Number(n)]).map(i=>Number.isFinite(i)?i:0);return e.length>=3?[e[0]||0,e[1]||0,e[1]||0,e[2]||0]:e.length===1?[e[0]||0,e[0]||0,e[0]||0,e[0]||0]:[0,0,0,0]};function to(n,t={}){var F,A,k,v,S,T;if(n.byteLength<8)return null;const e=new DataView(n.buffer,n.byteOffset,n.byteLength),i=e.getUint16(0,!1),r=i===18761;if(!r&&i!==19789||e.getUint16(2,r)!==42)return null;const s=e.getUint32(4,r);let o;try{o=Qr(e,r,s)}catch{return null}const a=new Set([s]);for(const L of Ut(e,r,o,330))L>0&&a.add(L);let l=null;for(const L of a)try{const E=L===s?o:Qr(e,r,L),V=Ut(e,r,E,256)[0]||0,U=Ut(e,r,E,257)[0]||0,N=Ut(e,r,E,258),D=Ut(e,r,E,277)[0]||N.length||1,B=N.length===1?new Array(D).fill(N[0]):N,z=Ut(e,r,E,259)[0]||1,q=Ut(e,r,E,262)[0]||0,O=Ut(e,r,E,284)[0]||1,W=Ut(e,r,E,339),j=W.length===0?new Array(D).fill(1):W.length===1?new Array(D).fill(W[0]):W,K=Ut(e,r,E,273),et=Ut(e,r,E,279),tt=Ut(e,r,E,278)[0]||U,ot=j.slice(0,D).every(R=>R===0||R===1);if(!(q===34892&&z===1&&O===1&&V>0&&U>0&&D>=3&&B.length>=D&&B.slice(0,D).every(R=>R===16)&&ot&&K.length>0&&K.length===et.length))continue;(!l||V*U>l.width*l.height)&&(l={entries:E,width:V,height:U,samplesPerPixel:D,rowsPerStrip:tt,stripOffsets:K,stripByteCounts:et,bitsPerSample:B,photometric:q})}catch{continue}if(!l)return null;const{entries:u,width:h,height:c,samplesPerPixel:y,rowsPerStrip:f,stripOffsets:p,stripByteCounts:d,bitsPerSample:g,photometric:x}=l,m=h*c;if(!Number.isSafeInteger(m)||m<=0)return null;const b=new Uint16Array(m*3);for(let L=0;L<p.length;L++){const E=L*f;if(E>=c)break;const V=Math.min(f,c-E),U=V*h*y*2,N=p[L],D=d[L];if(N<0||D<U||N+U>n.byteLength)throw new Error("Invalid LinearRaw DNG strip bounds");const B=E*h*3;if(y===3&&r&&!(n.byteOffset+N&1)){b.set(new Uint16Array(n.buffer,n.byteOffset+N,V*h*3),B);continue}const z=new DataView(n.buffer,n.byteOffset+N,U);let q=0,O=B;for(let W=0;W<V*h;W++)b[O++]=z.getUint16(q,r),b[O++]=z.getUint16(q+2,r),b[O++]=z.getUint16(q+4,r),q+=y*2}const _=Ut(e,r,u,50714),P=((A=(F=t==null?void 0:t.color_data)==null?void 0:F.dng_levels)==null?void 0:A.dng_cblack)||((k=t==null?void 0:t.color_data)==null?void 0:k.cblack_rawpy_style)||(t==null?void 0:t.black_level_per_channel)||(t==null?void 0:t.cblack),M=Za(_.length?_:P),w=Ut(e,r,u,50717),C=As(w,(S=(v=t==null?void 0:t.color_data)==null?void 0:v.dng_levels)==null?void 0:S.dng_whitelevel,(T=t==null?void 0:t.color_data)==null?void 0:T.maximum,t==null?void 0:t.white_level)||65535;return{data:b,width:h,height:c,bayerPattern:"",blackLevels:M,whiteLevel:C,metadata:{...t,format:"DNG_LINEAR_RAW_RGB",description:"Uncompressed DNG LinearRaw RGB",linearRawDngDecoder:!0,bitsPerSample:g.slice(0,y),samplesPerPixel:y,photometric:x},isThreePlane:!0,threePlaneTransfer:"linear"}}async function eo(n){const t=to(n);if(!t)return null;let e={};const i=await $i(),r=new i;try{await r.open(n,{}),e=await r.metadata(!0)}catch(s){console.warn("LinearRaw DNG metadata enrichment failed",s)}finally{try{r.delete?r.delete():r.close()}catch{}}try{const s=await Fs.parse(n.buffer);s&&(e={...e,...s})}catch(s){console.warn("exifr parsing failed for LinearRaw DNG",s)}return t.metadata={...e,...t.metadata},t}const no=async n=>{var r,s,o,a,l,u,h,c,y,f,p;const t=await eo(n);if(t)return t;const e=await $i(),i=new e;try{if(await i.open(n,{}),typeof i.getRawImage!="function")throw new Error("WASM mismatch");let d={};try{d=await i.metadata(!0)}catch(S){console.warn("Metadata error before raw extraction",S)}const g=i.getRawImage(),x=g.data instanceof Uint16Array?g.data:new Uint16Array(g.data);let m={...d};try{const S=await Fs.parse(n.buffer);S&&(m={...m,...S})}catch(S){console.warn("exifr parsing failed for RAW buffer",S)}const b=d.filters||((r=d.idata)==null?void 0:r.filters)||0,_=d.colors||((s=d.idata)==null?void 0:s.colors)||0,P=b===0&&_===3,M=b===9;let w=[0,0,0,0],C=null,F;if(i.getBlackLevels)try{const S=i.getBlackLevels();C=S,F=Tr(S)||Tr((o=d.color_data)==null?void 0:o.black_level_model)||void 0;const T=Qe(F==null?void 0:F.siteBaseLevels);T&&(w=T)}catch(S){console.warn("getBlackLevels binding failed",S)}if(!F){const S=Number((C==null?void 0:C.black)??((a=d.color_data)==null?void 0:a.black)??0)||0,T=Qe((C==null?void 0:C.cblack)||((l=d.color_data)==null?void 0:l.cblack_rawpy_style)||d.black_level_per_channel||d.cblack||((u=d.color)==null?void 0:u.cblack));if(T){const L=String(g.bayerPattern||d.cfa_pattern||"RGGB").toUpperCase();let E=0;w=[0,1,2,3].map(V=>{const U=L[V],N=U==="R"?0:U==="B"?2:++E===1?1:3;return Math.max(0,S+T[N])})}else w=[S,S,S,S]}const A=Number.isFinite(Number(g.bits))&&Number(g.bits)>0?Math.pow(2,Number(g.bits))-1:null,k=As(d.white_level,(c=(h=d.color_data)==null?void 0:h.dng_levels)==null?void 0:c.dng_whitelevel,(y=d.color_data)==null?void 0:y.maximum,C==null?void 0:C.maximum,(f=d.color_data)==null?void 0:f.fmaximum,(p=d.color_data)==null?void 0:p.data_maximum,A)||16383,v={data:x,width:g.width,height:g.height,bayerPattern:g.bayerPattern||"",blackLevels:w,blackLevelModel:F,whiteLevel:k,metadata:m,isThreePlane:P,threePlaneTransfer:P?"linear":void 0,isXTrans:M};return Ca(v,Pa(v,m)),v}finally{i.delete?i.delete():i.close()}};async function io(n){if(ga(n)){const i=await wa(n);if(!i)throw new Error("Sony cRAW HQ decoder did not return image data.");return i.rawImageData}const e=new Uint8Array(n);return no(e)}function ro(n,t,e){const i=Ts(n,e),r=Math.floor(t.x),s=Math.floor(t.y),o=Math.floor(t.w),a=Math.floor(t.h),l=new Uint16Array(o*a);for(let u=0;u<a;u++){const h=s+u,c=u*o;for(let y=0;y<o;y++)l[c+y]=i(r+y,h)}return{data:l,width:o,height:a}}function so(n,t,e){const i=ao(n,t,e);if(!i)return null;const r=n.width,s=n.height,o=new Uint16Array(r*s);for(let a=0;a<s;a++){const l=a*r;for(let u=0;u<r;u++)o[l+u]=i(u,a)}return{kind:"u16-mono",data:o,width:r,height:s}}function ao(n,t,e){return t.renderMode==="advanced-zero-dep"&&t.advancedZeroDep?Ts(n,t,e):t.renderMode==="zero-dependency"?oo(n,t,e):null}function Ts(n,t,e){if(!t.advancedZeroDep)throw new Error("Unmixing settings not found in DisplaySettings.");const{bg:i,fg:r}=t.advancedZeroDep,s=Is(e,t.advancedZeroDep.bl),{data:o,width:a,whiteLevel:l}=n,u=i.map((p,d)=>Math.max(0,p-s[d])),h=r.map((p,d)=>Math.max(0,p-s[d])),c=(u[1]+u[3])/2,y=(h[1]+h[3])/2,f=Math.pow(2,t.exposure);return(p,d)=>{if(p<0||d<0||p>=a||d>=n.height)return 0;const g=d%2,x=p%2;let m=0;!g&&!x?m=0:!g&&x?m=1:g&&!x?m=3:m=2;const b=o[d*a+p],_=s[m],P=u[m],M=h[m],w=Math.max(b-_,0),C=M-P||1e-9,F=(w-P)/C;let A;return F<0?A=w*(c/Math.max(P,1e-9)):F>1?A=w*(y/Math.max(M,1e-9)):A=(1-F)*c+F*y,A*=f,Math.max(0,Math.min(65535,Math.round(A)))}}function oo(n,t,e){const{data:i,width:r,height:s,whiteLevel:o}=n,a=Is(e,t.blackLevel||[0,0,0,0]),l=co(n.bayerPattern),u=t.wbGains?t.wbGains[0]:1,h=t.wbGains?t.wbGains[1]:1,c=Math.pow(2,t.exposure||0);return(y,f)=>{if(y<0||f<0||y>=r||f>=s)return 0;const p=uo(y,f),d=ho(l,y,f),g=i[f*r+y],x=a[p];let m=(g-x)/Math.max(1,o-x);return m=Math.max(0,Math.min(1,m)),m*=c,d==="R"?m*=u:d==="B"&&(m*=h),lo(m)}}function Is(n,t){if(typeof n=="number"&&Number.isFinite(n)){const e=Math.max(0,n);return[e,e,e,e]}return Array.isArray(n)&&n.length===4?[Math.max(0,n[0]??0),Math.max(0,n[1]??0),Math.max(0,n[2]??0),Math.max(0,n[3]??0)]:[Math.max(0,t[0]??0),Math.max(0,t[1]??0),Math.max(0,t[2]??0),Math.max(0,t[3]??0)]}function lo(n){return Math.max(0,Math.min(65535,Math.round(Math.max(0,Math.min(1,n))*65535)))}function co(n){const t=Ke(n);if(!t)throw new Error("Cannot process RAW mosaic without a valid Bayer CFA pattern.");return t}function uo(n,t){return ps(n,t)}function ho(n,t,e){return ka(n,t,e)}const ct=(n,t=0)=>({real:n,imag:t}),mn=(n,t)=>({real:n.real+t.real,imag:n.imag+t.imag}),dn=(n,t)=>({real:n.real-t.real,imag:n.imag-t.imag}),Pt=(n,t)=>({real:n.real*t.real-n.imag*t.imag,imag:n.real*t.imag+n.imag*t.real}),Ae=(n,t)=>{const e=t.real*t.real+t.imag*t.imag;return e===0?ct(0):{real:(n.real*t.real+n.imag*t.imag)/e,imag:(n.imag*t.real-n.real*t.imag)/e}},qe=n=>Math.hypot(n.real,n.imag),Ns=n=>{const t=qe(n);if(t===0)return ct(0);const e=Math.sqrt(t),i=Math.atan2(n.imag,n.real);return ct(e*Math.cos(i/2),e*Math.sin(i/2))};function fo(n,t){const e=n.length-1;if(e<0)return{p:ct(0),dp:ct(0),d2p:ct(0)};let i=ct(n[e].real,n[e].imag),r=ct(0),s=ct(0);for(let o=e-1;o>=0;o--)s=mn(Pt(r,ct(2)),Pt(t,s)),r=mn(i,Pt(t,r)),i=mn(ct(n[o].real,n[o].imag),Pt(t,i));return{p:i,dp:r,d2p:s}}function qi(n,t,e=80){const r=n.length-1;if(r<=0)return{root:t,iterations:0};if(r===1)return{root:Ae(Pt(n[0],ct(-1)),n[1]),iterations:0};let s=ct(t.real,t.imag);for(let o=0;o<e;o++){const{p:a,dp:l,d2p:u}=fo(n,s);if(qe(a)<1e-14)return{root:s,iterations:o};const h=Ae(l,a),c=Pt(h,h),y=dn(c,Ae(u,a)),f=ct(r),p=ct(r-1),d=dn(Pt(f,y),Pt(h,h)),g=Ns(Pt(p,d)),x=mn(h,g),m=dn(h,g),b=qe(x)>qe(m)?x:m;if(qe(b)<1e-14)return{root:s,iterations:o};const _=Ae(f,b),P=dn(s,_);if(qe(_)<1e-14*qe(P))return{root:P,iterations:o+1};s=P}return{root:s,iterations:e}}function po(n,t){const e=n.length-1;if(e<=0)return[ct(0)];if(e===1)return[ct(n[0].real,n[0].imag)];const i=new Array(e);i[e-1]=ct(n[e].real,n[e].imag);for(let r=e-2;r>=0;r--){const s=ct(n[r+1].real,n[r+1].imag),o=i[r+1];i[r]=mn(s,Pt(t,o))}return i}function mo(n){const t=n.length-1;if(t<=0)return[];if(t===1)return[Ae(Pt(n[0],ct(-1)),n[1])];const e=[];let i=n.map(s=>ct(s.real,s.imag)),r=t*5;for(;i.length>2&&r-- >0;){const s=ct(.3+Math.random()*.7,.3+Math.random()*.7),{root:o}=qi(i,s,100),a=qi(n,o,20);e.push(a.root);const l=po(i,o);if(l.length>=i.length){console.warn("polyDeflate did not reduce degree, breaking");break}i=l}if(i.length===2)e.push(Ae(Pt(i[0],ct(-1)),i[1]));else if(i.length===3){const s=i[2],o=i[1],a=i[0],l=dn(Pt(o,o),Pt(Pt(ct(4),s),a)),u=Ns(l),h=Pt(ct(2),s),c=Ae(dn(Pt(ct(-1),o),u),h),y=Ae(mn(Pt(ct(-1),o),u),h);e.push(c,y)}return e}function go(n,t,e){const i=[ct(n),ct(-1),ct(n*t),ct(0),ct(n*e)],r=mo(i);if(r.length===0)return console.warn("laguerreSmallestPositiveRoot: no roots found"),n;let s=1/0,o=!1;for(const l of r)Math.abs(l.imag)<1e-10&&l.real>0&&l.real<s&&(s=l.real,o=!0);return o?qi(i,ct(s,0),20).root.real:(console.warn("laguerreSmallestPositiveRoot: no positive real root found"),n)}function yo(n,t,e){if(Math.abs(t)<1e-10&&Math.abs(e)<1e-10)return n;if(n<1e-10)return 0;if(Math.abs(e)<1e-10){const i=-1/(t*n),r=1/t,s=i*i-4*r;if(s<0)return n;const o=Math.sqrt(s),a=-.5*(i+Math.sign(i)*o),l=a,u=r/a;return l>0&&u>0?Math.min(l,u):l>0?l:u>0?u:n}try{return go(n,t,e)}catch(i){return console.error("laguerreSmallestPositiveRoot failed:",i),n}}function xo(n,t,e,i,r){const s=n.x-t.x,o=n.y-t.y,a=Math.hypot(s,o)/Math.max(1e-12,e);if(a<1e-12)return{x:n.x,y:n.y};const l=a*a,u=1+(i+r*l)*l;return{x:s/u+t.x,y:o/u+t.y}}function bo(n,t,e,i,r){const s=n.x-t.x,o=n.y-t.y,a=Math.hypot(s,o)/Math.max(1e-12,e);if(a<1e-12)return{x:t.x,y:t.y};const u=yo(a,i,r)/a;return{x:t.x+s*u,y:t.y+o*u}}function Tt(n,t){const e=xo(n,{x:t.principalX,y:t.principalY},t.radiusNorm,t.k1,t.k2);return{x:e.x+(t.correctedOffsetX??0),y:e.y+(t.correctedOffsetY??0)}}function we(n,t){const e={x:n.x-(t.correctedOffsetX??0),y:n.y-(t.correctedOffsetY??0)};return bo(e,{x:t.principalX,y:t.principalY},t.radiusNorm,t.k1,t.k2)}const _o=`
attribute vec2 a_position;

void main() {
  gl_Position = vec4(a_position, 0.0, 1.0);
}
`,wo=`
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
`,Mo=`
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
`;class So{constructor(){yt(this,"canvas",null);yt(this,"gl",null);yt(this,"blurProgram",null);yt(this,"sobelProgram",null);yt(this,"positionBuffer",null);yt(this,"blurUniforms",null);yt(this,"sobelUniforms",null);yt(this,"resources",null);yt(this,"initialized",!1);yt(this,"unavailable",!1);yt(this,"maxTextureSize",0)}compute(t,e,i){if(!this.initialized&&!this.init())return null;const r=this.gl,s=this.blurProgram,o=this.sobelProgram,a=this.blurUniforms,l=this.sobelUniforms;if(!r||!s||!o||!a||!l||!this.positionBuffer||!this.canvas||e<=2||i<=2||e>this.maxTextureSize||i>this.maxTextureSize)return null;const u=this.ensureResources(e,i);if(!u)return null;this.canvas.width=e,this.canvas.height=i,r.viewport(0,0,e,i),r.disable(r.BLEND),r.pixelStorei(r.UNPACK_ALIGNMENT,1),r.pixelStorei(r.PACK_ALIGNMENT,1),r.bindBuffer(r.ARRAY_BUFFER,this.positionBuffer),r.activeTexture(r.TEXTURE0),r.bindTexture(r.TEXTURE_2D,u.sourceTexture),r.texImage2D(r.TEXTURE_2D,0,r.LUMINANCE,e,i,0,r.LUMINANCE,r.UNSIGNED_BYTE,t),r.useProgram(s),r.enableVertexAttribArray(0),r.vertexAttribPointer(0,2,r.FLOAT,!1,0,0),r.uniform2f(a.size,e,i),r.uniform1i(a.source,0),r.bindFramebuffer(r.FRAMEBUFFER,u.blurFramebuffer),r.bindTexture(r.TEXTURE_2D,u.sourceTexture),r.drawArrays(r.TRIANGLES,0,6);const h=new Uint8Array(e*i*4);r.readPixels(0,0,e,i,r.RGBA,r.UNSIGNED_BYTE,h),r.useProgram(o),r.uniform2f(l.size,e,i),r.uniform1i(l.blurred,0),r.bindFramebuffer(r.FRAMEBUFFER,u.sobelFramebuffer),r.bindTexture(r.TEXTURE_2D,u.blurTexture),r.drawArrays(r.TRIANGLES,0,6);const c=new Uint8Array(e*i*4);r.readPixels(0,0,e,i,r.RGBA,r.UNSIGNED_BYTE,c),r.disableVertexAttribArray(0),r.bindFramebuffer(r.FRAMEBUFFER,null),r.bindBuffer(r.ARRAY_BUFFER,null),r.bindTexture(r.TEXTURE_2D,null);const y=new Uint8Array(e*i);for(let g=0,x=0;g<y.length;g++,x+=4)y[g]=h[x];const f=new Float32Array(e*i),p=new Float32Array(e*i),d=new Float32Array(e*i);for(let g=0,x=0;g<f.length;g++,x+=4){const m=(c[x]|c[x+1]<<8)-32768,b=(c[x+2]|c[x+3]<<8)-32768;f[g]=m,p[g]=b,d[g]=Math.sqrt(m*m+b*b)}return{blurredGray:y,gx:f,gy:p,magnitude:d}}init(){if(this.initialized&&this.gl&&this.blurProgram&&this.sobelProgram)return!0;if(this.unavailable)return!1;const t=this.createCanvas();if(!t)return this.unavailable=!0,!1;const e=t.getContext("webgl",{alpha:!1,antialias:!1,depth:!1,stencil:!1,premultipliedAlpha:!1,preserveDrawingBuffer:!1});if(!e)return this.unavailable=!0,!1;const i=this.compileShader(e,e.VERTEX_SHADER,_o),r=this.compileShader(e,e.FRAGMENT_SHADER,wo),s=this.compileShader(e,e.FRAGMENT_SHADER,Mo);if(!i||!r||!s)return i&&e.deleteShader(i),r&&e.deleteShader(r),s&&e.deleteShader(s),this.unavailable=!0,!1;const o=this.createProgram(e,i,r),a=this.createProgram(e,i,s);if(e.deleteShader(i),e.deleteShader(r),e.deleteShader(s),!o||!a)return o&&e.deleteProgram(o),a&&e.deleteProgram(a),this.unavailable=!0,!1;const l=e.createBuffer();return l?(e.bindBuffer(e.ARRAY_BUFFER,l),e.bufferData(e.ARRAY_BUFFER,new Float32Array([-1,-1,1,-1,-1,1,-1,1,1,-1,1,1]),e.STATIC_DRAW),e.bindBuffer(e.ARRAY_BUFFER,null),this.canvas=t,this.gl=e,this.blurProgram=o,this.sobelProgram=a,this.positionBuffer=l,this.blurUniforms={source:e.getUniformLocation(o,"u_source"),size:e.getUniformLocation(o,"u_size")},this.sobelUniforms={blurred:e.getUniformLocation(a,"u_blurred"),size:e.getUniformLocation(a,"u_size")},this.maxTextureSize=Number(e.getParameter(e.MAX_TEXTURE_SIZE)||0),this.initialized=!0,!0):(e.deleteProgram(o),e.deleteProgram(a),this.unavailable=!0,!1)}createCanvas(){return typeof OffscreenCanvas<"u"?new OffscreenCanvas(1,1):typeof document<"u"?document.createElement("canvas"):null}ensureResources(t,e){const i=this.gl;if(!i)return null;if(this.resources&&this.resources.width===t&&this.resources.height===e)return this.resources;this.disposeResources();const r=this.createTexture(i.LUMINANCE,t,e,i.LUMINANCE,i.UNSIGNED_BYTE,null),s=this.createTexture(i.RGBA,t,e,i.RGBA,i.UNSIGNED_BYTE,null),o=this.createTexture(i.RGBA,t,e,i.RGBA,i.UNSIGNED_BYTE,null),a=this.createFramebuffer(s),l=this.createFramebuffer(o);return!r||!s||!o||!a||!l?(r&&i.deleteTexture(r),s&&i.deleteTexture(s),o&&i.deleteTexture(o),a&&i.deleteFramebuffer(a),l&&i.deleteFramebuffer(l),null):(this.resources={width:t,height:e,sourceTexture:r,blurTexture:s,sobelTexture:o,blurFramebuffer:a,sobelFramebuffer:l},this.resources)}createTexture(t,e,i,r,s,o){const a=this.gl;if(!a)return null;const l=a.createTexture();return l?(a.bindTexture(a.TEXTURE_2D,l),a.texParameteri(a.TEXTURE_2D,a.TEXTURE_MIN_FILTER,a.NEAREST),a.texParameteri(a.TEXTURE_2D,a.TEXTURE_MAG_FILTER,a.NEAREST),a.texParameteri(a.TEXTURE_2D,a.TEXTURE_WRAP_S,a.CLAMP_TO_EDGE),a.texParameteri(a.TEXTURE_2D,a.TEXTURE_WRAP_T,a.CLAMP_TO_EDGE),a.texImage2D(a.TEXTURE_2D,0,t,e,i,0,r,s,o),a.bindTexture(a.TEXTURE_2D,null),l):null}createFramebuffer(t){const e=this.gl;if(!e||!t)return null;const i=e.createFramebuffer();if(!i)return null;e.bindFramebuffer(e.FRAMEBUFFER,i),e.framebufferTexture2D(e.FRAMEBUFFER,e.COLOR_ATTACHMENT0,e.TEXTURE_2D,t,0);const r=e.checkFramebufferStatus(e.FRAMEBUFFER);return e.bindFramebuffer(e.FRAMEBUFFER,null),r!==e.FRAMEBUFFER_COMPLETE?(e.deleteFramebuffer(i),null):i}compileShader(t,e,i){const r=t.createShader(e);return r?(t.shaderSource(r,i),t.compileShader(r),t.getShaderParameter(r,t.COMPILE_STATUS)?r:(console.error("[SFR Auto Detect WebGL] shader compile failed",t.getShaderInfoLog(r)),t.deleteShader(r),null)):null}createProgram(t,e,i){const r=t.createProgram();return r?(t.attachShader(r,e),t.attachShader(r,i),t.bindAttribLocation(r,0,"a_position"),t.linkProgram(r),t.getProgramParameter(r,t.LINK_STATUS)?r:(console.error("[SFR Auto Detect WebGL] program link failed",t.getProgramInfoLog(r)),t.deleteProgram(r),null)):null}disposeResources(){const t=this.gl,e=this.resources;if(!t||!e){this.resources=null;return}t.deleteTexture(e.sourceTexture),t.deleteTexture(e.blurTexture),t.deleteTexture(e.sobelTexture),t.deleteFramebuffer(e.blurFramebuffer),t.deleteFramebuffer(e.sobelFramebuffer),this.resources=null}}const vo=new So,Po=12,Co=12,ko=501,Fo=2,Ao=[.3,.4,.55,.75,1,1.35],To=320,Io=.75,No=.002,Kr=32,Bt=1e-12;function Ht(n,t,e){return Math.max(t,Math.min(e,n))}function fn(n){if(n.length===0)return NaN;const t=[...n].sort((i,r)=>i-r),e=t.length>>1;return t.length&1?t[e]:.5*(t[e-1]+t[e])}function $r(n,t){if(n.length===0)return NaN;const e=Ht(t,0,1)*(n.length-1),i=Math.floor(e),r=Math.ceil(e),s=e-i;return n[i]+s*(n[r]-n[i])}function ar(n){if(n.length===0)return 0;const t=fn(n),e=n.map(i=>Math.abs(i-t));return 1.4826*fn(e)}function Ro(n,t){const e=n.map(a=>(a.x%1+1)%1).sort((a,l)=>a-l),i=[];for(const a of e)(i.length===0||a-i[i.length-1]>1e-5)&&i.push(a);let r=1;if(i.length>1){r=1-i[i.length-1]+i[0];for(let a=1;a<i.length;a++)r=Math.max(r,i[a]-i[a-1])}const s=Ht(1-r,0,1);let o=t;return i.length<8||s<.75?o=Math.min(o,.5):s<.9&&(o=Math.min(o,1)),n.length<96&&(o=Math.min(o,.75)),{uniquePhaseCount:i.length,maximumPhaseGapPx:r,phaseCoverage:s,reliableFrequency:o}}function Lo(n,t){let e=0,i=n.length;for(;e<i;){const r=e+i>>1;n[r].x<t?e=r+1:i=r}return e}function Eo(n,t){let e=0,i=n.length;for(;e<i;){const r=e+i>>1;n[r].x<=t?e=r+1:i=r}return e}function or(n,t){const e=t.length,i=n.map((r,s)=>[...r,t[s]]);for(let r=0;r<e;r++){let s=r;for(let a=r+1;a<e;a++)Math.abs(i[a][r])>Math.abs(i[s][r])&&(s=a);if(Math.abs(i[s][r])<=1e-14)return null;s!==r&&([i[s],i[r]]=[i[r],i[s]]);const o=i[r][r];for(let a=r;a<=e;a++)i[r][a]/=o;for(let a=0;a<e;a++){if(a===r)continue;const l=i[a][r];if(l!==0)for(let u=r;u<=e;u++)i[a][u]-=l*i[r][u]}}return i.map(r=>r[e])}function Rs(n,t,e,i,r,s=-1){const o=Lo(n,t-e),a=Eo(n,t+e),l=i+1,u=Array.from({length:l},()=>new Array(l).fill(0)),h=new Array(l).fill(0);let c=0,y=0;for(let d=o;d<a;d++){const g=n[d];if(s>=0&&g.fold===s)continue;const x=(g.x-t)/e,m=Math.abs(x);if(m>=1)continue;const _=Math.pow(1-m*m*m,3)*(r?r[d]:1);if(!(_>0)||!Number.isFinite(_))continue;const P=new Array(l).fill(1);for(let M=1;M<l;M++)P[M]=P[M-1]*x;for(let M=0;M<l;M++){h[M]+=_*P[M]*g.y;for(let w=0;w<l;w++)u[M][w]+=_*P[M]*P[w]}c++,y+=_}if(c<Math.max(8,l*2)||y<=Bt)return null;const f=Math.max(Bt,u[0][0]*1e-11);for(let d=0;d<l;d++)u[d][d]+=f*Math.pow(4,d);const p=or(u,h);return!p||p.some(d=>!Number.isFinite(d))?null:{value:p[0],derivative:p[1]/e,effectiveCount:c,weightSum:y}}function Uo(n,t){const e=t*.65,i=n.filter(o=>Math.abs(o.x)>=e);if(i.length<12)return 0;let r=new Array(i.length).fill(1),s=null;for(let o=0;o<3;o++){const a=Array.from({length:3},()=>new Array(3).fill(0)),l=new Array(3).fill(0);for(let c=0;c<i.length;c++){const y=i[c],f=[1,y.x>=0?1:0,y.x/t],p=r[c];for(let d=0;d<3;d++){l[d]+=p*f[d]*y.y;for(let g=0;g<3;g++)a[d][g]+=p*f[d]*f[g]}}if(s=or(a,l),!s)return 0;const u=i.map(c=>{const y=s[0]+s[1]*(c.x>=0?1:0)+s[2]*c.x/t;return c.y-y}),h=Math.max(Bt,ar(u));r=u.map(c=>{const y=Math.abs(c)/(1.345*h);return y<=1?1:1/y})}return s?s[2]/t:0}function Do(n,t,e,i){const r=[];for(let m=0;m<n.length;m++){const b=Number(n[m]),_=Number(t[m]);Number.isFinite(b)&&Number.isFinite(_)&&r.push({x:b,y:_,fold:m%5})}if(r.sort((m,b)=>m.x-b.x),r.length<Kr)return null;const s=Math.min(Math.abs(r[0].x),Math.abs(r[r.length-1].x)),o=Ht(Number.isFinite(e)?e:Math.min(Co,s*.96),3,Math.max(3,s*.98));let a=r.filter(m=>Math.abs(m.x)<=o);if(a.length<Kr)return null;const l=Ht(Math.floor(a.length*.08),8,64),u=fn(a.slice(0,l).map(m=>m.y)),h=fn(a.slice(-l).map(m=>m.y));if(!Number.isFinite(u)||!Number.isFinite(h)||Math.abs(h-u)<=Bt)return null;u>h&&(a=a.map(m=>({...m,x:-m.x})).sort((m,b)=>m.x-b.x));const c=i?Uo(a,o):0,y=a.map(m=>({...m,y:m.y-c*m.x})),f=fn(y.slice(0,l).map(m=>m.y)),d=fn(y.slice(-l).map(m=>m.y))-f;if(!Number.isFinite(d)||d<=Bt)return null;const g=y.map(m=>({...m,y:(m.y-f)/d})),x=[...g.slice(0,l).map(m=>m.y),...g.slice(-l).map(m=>m.y-1)];return{samples:g,halfWidth:o,normalization:{dark:f,contrast:d,removedSlopePerPixel:c,noiseSigmaNormalized:ar(x)}}}function Bo(n,t){const e=t.map(Number).filter(h=>Number.isFinite(h)&&h>=.2&&h<=2.5).sort((h,c)=>h-c);if(e.length===0)return .55;const i=n.filter(h=>Math.abs(h.x)<=4),r=i.length>=64?i:n,s=Math.min(320,r.length),o=[];for(let h=0;h<s;h++)o.push(r[Math.floor((h+.5)*r.length/s)]);const a=[];for(const h of e){const c=[];for(const y of o){const f=Rs(n,y.x,h,3,null,y.fold);f&&c.push(Math.abs(f.value-y.y))}if(c.length>=o.length*.8){const y=c.sort((d,g)=>d-g),f=$r(y,.5),p=$r(y,.8);a.push({bandwidth:h,score:f+.2*p})}}if(a.length===0)return e[e.length-1];a.sort((h,c)=>h.score-c.score);const l=a[0];return a.filter(h=>h.score<=l.score*1.025+1e-8).reduce((h,c)=>Math.max(h,c.bandwidth),l.bandwidth)}function Ls(n,t,e){if(n<=t[0])return e[0];if(n>=t[t.length-1])return e[e.length-1];let i=0,r=t.length-1;for(;r-i>1;){const o=i+r>>1;t[o]<=n?i=o:r=o}const s=(n-t[i])/(t[r]-t[i]);return e[i]+s*(e[r]-e[i])}function Jr(n,t,e,i,r,s=null){const o=1/e,a=Math.max(33,Math.round(2*t/o)+1),l=new Array(a);for(let d=0;d<a;d++)l[d]=-t+d*(2*t/(a-1));const u=s?n.map((d,g)=>({...d,y:s(d,g)})):n;let h=new Array(u.length).fill(1),c=[],y=[],f=[];const p=()=>{c=new Array(a),y=new Array(a),f=new Array(a);for(let d=0;d<a;d++){let g=i,x=null;for(let m=0;m<4&&!x;m++)x=Rs(u,l[d],g,3,h),g*=1.35;if(!x)return!1;c[d]=x.value,y[d]=x.derivative,f[d]=x.effectiveCount}return!0};if(!p())return null;if(!s)for(let d=0;d<r;d++){const g=u.map(m=>m.y-Ls(m.x,l,c)),x=Math.max(Bt,ar(g));if(h=g.map(m=>{const b=m/(4.685*x);if(Math.abs(b)>=1)return .001;const _=1-b*b;return Math.max(.001,_*_)}),!p())return null}return{x:l,esf:c,lsf:y,supportCounts:f,step:l[1]-l[0]}}function Oo(n){let t=1;for(;t<n;)t*=2;return t}function zo(n,t=null){const e=n.length;if(e<1||e&e-1)throw new Error("FFT length must be a power of two.");const i=Float64Array.from(n),r=t?Float64Array.from(t):new Float64Array(e);for(let s=1,o=0;s<e;s++){let a=e>>1;for(;o&a;)o^=a,a>>=1;o^=a,s<o&&([i[s],i[o]]=[i[o],i[s]],[r[s],r[o]]=[r[o],r[s]])}for(let s=2;s<=e;s*=2){const o=-2*Math.PI/s,a=Math.cos(o),l=Math.sin(o);for(let u=0;u<e;u+=s){let h=1,c=0;const y=s>>1;for(let f=0;f<y;f++){const p=u+f,d=p+y,g=i[d]*h-r[d]*c,x=i[d]*c+r[d]*h,m=i[p],b=r[p];i[p]=m+g,r[p]=b+x,i[d]=m-g,r[d]=b-x;const _=h*a-c*l;c=h*l+c*a,h=_}}}return{real:i,imag:r}}function Vo(n,t,e,i){const r=Math.abs(n-t)/Math.max(Bt,e);if(r>=1)return 0;if(r<=i)return 1;const s=(r-i)/(1-i);return .5*(1+Math.cos(Math.PI*s))}function Zr(n,t,e=null,i=null){const r=n.lsf.reduce((m,b,_)=>Math.abs(b)>Math.abs(n.lsf[m])?_:m,0),s=n.x[r],o=Number.isFinite(i)?i:s,a=Math.min(o-n.x[0],n.x[n.x.length-1]-o),l=Math.max(2,a*.98),u=n.lsf.map((m,b)=>m*Vo(n.x[b],o,l,.72)),h=e||Oo(Math.max(2048,u.length*8)),c=new Float64Array(h);c.set(u.slice(0,h));const y=zo(c),f=Math.max(Bt,Math.hypot(y.real[0],y.imag[0])),p=1/(h*n.step),d=Math.min(Math.floor(t/p)+1,h>>1),g=new Array(d),x=new Array(d);for(let m=0;m<d;m++)g[m]=m*p,x[m]=Math.hypot(y.real[m],y.imag[m])/f;return{frequencies:g,mtf:x,fftSize:h,center:o,peakPosition:s,windowHalfWidth:l}}function Go(n,t,e){for(let i=1;i<t.length;i++)if(t[i-1]>e&&t[i]<=e){const r=t[i]-t[i-1],s=Math.abs(r)<=Bt?1:(e-t[i-1])/r;return{frequency:n[i-1]+s*(n[i]-n[i-1]),index:i}}return null}function Xo(n,t,e,i){if(!e)return null;const r=n[1]-n[0],s=Math.max(8,Math.round(.035/Math.max(Bt,r))),o=Math.max(0,e.index-s),a=Math.min(t.length-1,e.index+s);if(a-o+1<9)return null;const l=e.frequency,u=Math.max(r,n[a]-n[o]),h=Array.from({length:3},()=>new Array(3).fill(0)),c=new Array(3).fill(0);for(let m=o;m<=a;m++){if(!(t[m]>0)||!Number.isFinite(t[m]))continue;const b=(n[m]-l)/u,_=Math.exp(-4*b*b),P=[1,b,b*b],M=Math.log(t[m]);for(let w=0;w<3;w++){c[w]+=_*P[w]*M;for(let C=0;C<3;C++)h[w][C]+=_*P[w]*P[C]}}const y=or(h,c);if(!y||y[1]>=0)return null;const f=y[0]-Math.log(i),p=y[1],d=y[2],g=[];if(Math.abs(d)<=1e-10)Math.abs(p)>Bt&&g.push(-f/p);else{const m=p*p-4*d*f;if(m>=0){const b=Math.sqrt(m);g.push((-p-b)/(2*d),(-p+b)/(2*d))}}return g.map(m=>l+m*u).filter(m=>Number.isFinite(m)&&m>=n[o]&&m<=n[a]).sort((m,b)=>Math.abs(m-l)-Math.abs(b-l))[0]??null}function zi(n,t,e){return n.map(i=>Ls(i,t,e))}function Yo(n,t,e){for(let i=1;i<t.length;i++)if(t[i-1]>=e&&t[i]<e){const r=t[i]-t[i-1],s=Math.abs(r)<=Bt?0:(e-t[i-1])/r;return n[i-1]+s*(n[i]-n[i-1])}return n[n.length-1]}function Es(n,t,e={}){if(!n||!t||n.length!==t.length)throw new TypeError("distancesPx and values must have the same length.");const i=Ht(Math.round(Number(e.oversampling)||Po),4,16),r=Ht(Number(e.maxFrequencyCyclesPerPixel)||Fo,.5,i*.45),s=Ht(Math.round(Number(e.outputPoints)||ko),65,4001),o=e.detrend!==!1,a=Number(e.halfWidthPx),l=Do(n,t,a,o);if(!l)return null;const u=l.samples.length<To,h=Number(e.robustIterations),c=Number.isFinite(h)?Ht(Math.round(h),0,4):l.normalization.noiseSigmaNormalized>=No?u?4:2:0,y=Number(e.bandwidthPx),f=Number.isFinite(y)?Ht(y,.2,2.5):u?Io:Bo(l.samples,e.bandwidthCandidatesPx||Ao),p=Jr(l.samples,l.halfWidth,i,f,c);if(!p)return null;const d=Jr(l.samples,l.halfWidth,i,f,0,z=>z.x>=0?1:0);if(!d)return null;const g=Zr(p,r,null,0),x=Zr(d,r,g.fftSize,0),m=Ht(Number(e.correctionFloor)||.18,.05,.8),b=Ht(Number(e.correctionMaxGain)||3,1,8),_=g.mtf.map((z,q)=>{const O=Math.max(m,x.mtf[q]||m);return z*Math.min(b,1/O)}),P=Ro(l.samples,r),M=Yo(x.frequencies,x.mtf,1/b),w=Math.min(P.reliableFrequency,M),C=Ht(Number(e.mtfLevel)||.5,.05,.95),F=Go(g.frequencies,_,C),A=Xo(g.frequencies,_,F,C),k=new Array(s);for(let z=0;z<s;z++)k[z]=z*r/(s-1);const v=zi(k,g.frequencies,g.mtf),S=zi(k,g.frequencies,_),T=zi(k,x.frequencies,x.mtf),L=Number(e.pixelPitchUm),E=Number.isFinite(L)&&L>0,V=E?1e3/L:1,U=k.map(z=>z*V),N=A??(F==null?void 0:F.frequency)??null,D=Math.min(...T.filter(Number.isFinite)),B=[];return Math.abs(g.peakPosition)>1&&B.push("LSF peak is more than one pixel from the supplied edge origin."),D<m&&B.push("High-frequency estimator correction reached its configured floor."),P.phaseCoverage<.75?B.push("Projected samples do not cover a complete subpixel phase cycle; frequencies above 0.5 cycles/pixel are not reliable."):P.phaseCoverage<.9&&B.push("Projected samples have incomplete subpixel phase coverage; high-frequency results have reduced confidence."),M<r-Bt&&B.push("Estimator correction reached its configured gain limit before the requested maximum frequency."),N!==null&&N>w+Bt&&B.push("The requested MTF crossing is above the reliable frequency limit for this sample geometry."),F||B.push(`MTF did not cross ${C.toFixed(3)} inside the requested frequency range.`),{esf:p.esf,lsf:p.lsf,profilePositionsPx:p.x,frequencies:U,mtf:S,mtfRaw:v,estimatorResponse:T,mtf50:N===null?null:N*V,mtf50Linear:F?F.frequency*V:null,mtf50Method:A===null?"linear":"local-log-quadratic",frequencyUnit:E?"cycles/mm":"cycles/pixel",diagnostics:{sampleCount:l.samples.length,oversampling:i,halfWidthPx:l.halfWidth,bandwidthPx:f,smallWindowMode:u,robustIterations:c,profileStepPx:p.step,fftSize:g.fftSize,lsfPeakPositionPx:g.peakPosition,windowHalfWidthPx:g.windowHalfWidth,correctionFloor:m,correctionMaxGain:b,samplingUniquePhaseCount:P.uniquePhaseCount,samplingMaximumPhaseGapPx:P.maximumPhaseGapPx,samplingPhaseCoverage:P.phaseCoverage,estimatorReliableFrequencyCyclesPerPixel:M,maxReliableFrequencyCyclesPerPixel:w,normalization:l.normalization,minLocalSupport:Math.min(...p.supportCounts),warnings:B}}}const cn=(n,t)=>{const e=Ke(n);if(!e)throw new Error(`Bayer CFA pattern is unresolved for ${t}.`);return e},Qi={gradientPercentiles:[.82,.88,.92,.95,.98,.995],downsampleMaxSide:1600,minComponentAreaRatio:15e-6,maxComponentAreaRatio:.35,minComponentAreaPx:20,minEdgePoints:24,extentQuantileLow:.02,extentQuantileHigh:.98,cornerTrimRatio:.18,minSpanPx:8,maxAspectRatio:2,bandScale:.16,bandMinPx:1.75,bandMaxPx:14,minPointContrast:6,minSidePoints:3,minCoverageRatio:.15,minCenterCoverageRatio:.2,filterBlockPurity:!0,innerPurityStdScale:1.5,outerMeanSpreadLimit:51,minAxisDot:.6,residualLimitFloor:.01,residualLimitScale:.25,minQuadArea:48,minSideLength:10,minOuterContrast:5,sampleHalfWidthRatio:.25};function Wo(n,t,e,i,r,s){const o=n.width,a=n.height,l=cn(n.bayerPattern,"corrected RAW SFR sampling"),u=[],h=[],c=s!=null&&s.correctedRect?Ft*2:Ft,y=Math.max(1,Math.min(r,c)),f=e.p2.x-e.p1.x,p=e.p2.y-e.p1.y,d=Math.hypot(f,p);if(!Number.isFinite(d)||d<=1e-6)return null;const g=f/d,x=p/d,m=-x,b=g,_={x:(e.p1.x+e.p2.x)*.5,y:(e.p1.y+e.p2.y)*.5},P=s!=null&&s.correctedRect?xt(s.correctedRect,o,a):xt(Nt(he(e,c*4+2)??[e.p1,e.p2],2),o,a);if(!P)return null;const M=(s==null?void 0:s.correctedScanlinesOverride)??(s!=null&&s.distortedRect?gr(xt(s.distortedRect,o,a)??s.distortedRect,t,o,a):Qs(P,e,Math.max(1,i),y*4+.5,o,a));if(!M||M.size===0)return null;const w=Ks(M,t,o,a);if(w.size===0)return null;const C=!yr(t);for(const[A,k]of w){if(A<0||A>=a)continue;const v=A*o;for(let S=k.start;S<=k.end;S++){if(S<0||S>=o||!_t(S,A,l,s==null?void 0:s.greenPhase))continue;const T={x:S,y:A},L=Tt(T,t);if(!Number.isFinite(L.x)||!Number.isFinite(L.y)||Math.round(L.x)<0||Math.round(L.x)>=o||Math.round(L.y)<0||Math.round(L.y)>=a)continue;const E=L.x-_.x,V=L.y-_.y,U=E*g+V*x;let N=E*m+V*b;if(C){const D=$s(U,g,x,_,T,t);if(!D)continue;const B=.5*(D.a+D.b),z=Xt(D.a,g,x,_,t),q=Xt(B,g,x,_,t),O=Xt(D.b,g,x,_,t),W=br({x:D.a,y:Math.hypot(z.x-T.x,z.y-T.y)},{x:B,y:Math.hypot(q.x-T.x,q.y-T.y)},{x:D.b,y:Math.hypot(O.x-T.x,O.y-T.y)});if(!Number.isFinite(W))continue;const j=xr(W,g,x,_,t),K=Math.hypot(j.x,j.y);if(!Number.isFinite(K)||K<=1e-9)continue;const et=j.x/K,ot=-(j.y/K),rt=et,R=Xt(W,g,x,_,t);N=(T.x-R.x)*ot+(T.y-R.y)*rt}!Number.isFinite(U)||Math.abs(U)>Math.max(1,i)||!Number.isFinite(N)||Math.abs(N)>y||(u.push(N),h.push(Math.max(0,n.data[v+S]-bn(s==null?void 0:s.blackLevel,S,A))))}}if(u.length<8)return null;const F=Math.abs(f)>=Math.abs(p)?1:2;return s!=null&&s.forceLegacyModel?wn(u,h,F,c):_n(u,h,F,c)}function lr(n,t){const e=n.length;let i=0,r=0,s=0,o=0;for(let l=0;l<e;l++)i+=n[l],r+=t[l],s+=n[l]*t[l],o+=n[l]*n[l];const a=e*o-i*i;return a===0?{slope:0,intercept:0}:{slope:(e*s-i*r)/a,intercept:(r*o-i*s)/a}}function jo(n,t){const e=n.length,i=new Array(e).fill(0),r=2*t;for(let s=0;s<e;s++){const o=s>0?n[s-1]:n[0],a=s<e-1?n[s+1]:n[e-1];i[s]=(a-o)/r}return i}const An=-1e7,si=13,ce=512,At=8,cr=1/At,Ft=28,Ho=[[0,0,0,0,0,-.085714285714286,.342857142857143,.485714285714286,.342857142857143,-.085714285714286,0,0,0,0,0],[0,0,0,0,-.095238095238095,.142857142857143,.285714285714286,.333333333333333,.285714285714286,.142857142857143,-.095238095238095,0,0,0,0],[0,0,0,-.090909090909091,.060606060606061,.168831168831169,.233766233766234,.255411255411255,.233766233766234,.168831168831169,.060606060606061,-.090909090909091,0,0,0],[0,0,-.083916083916084,.020979020979021,.102564102564103,.160839160839161,.195804195804196,.207459207459208,.195804195804196,.160839160839161,.102564102564103,.020979020979021,-.083916083916084,0,0],[0,-.076923076923077,0,.062937062937063,.111888111888112,.146853146853147,.167832167832168,.174825174825175,.167832167832168,.146853146853147,.111888111888112,.062937062937063,0,-.076923076923077,0],[-.070588235294118,-.011764705882353,.038009049773756,.078733031674208,.110407239819004,.133031674208145,.146606334841629,.151131221719457,.146606334841629,.133031674208145,.110407239819004,.078733031674208,.038009049773756,-.011764705882353,-.070588235294118]];function ai(n,t,e,i=1){const r=Math.max(1e-6,e*.5),s=Math.max(1e-6,t*i),o=Math.exp(-s*r),a=1-o;if(!Number.isFinite(a)||Math.abs(a)<=1e-9)return Math.abs(n)<=r?1:0;if(Math.abs(n)<r){const l=2-2*o*Math.cosh(s*n),u=2*Math.sinh(s*r)*a;return!Number.isFinite(u)||Math.abs(u)<=1e-9?0:l/u}return Math.exp(-s*Math.abs(n))/a}function qo(n,t,e,i,r,s){const o=n.length;if(o===0)return[];if(s<1)return n;const a=Math.min(s,32),l=new Array(o).fill(0);l[0]=n[0];for(let f=1;f<o;f++)l[f]=l[f-1]+n[f];const u=(f,p)=>{const d=Math.max(0,f),g=Math.min(o-1,p);return g<d?n[Math.max(0,Math.min(o-1,f))]??0:(l[g]-(d>0?l[d-1]:0))/(g-d+1)},h=a*2,c=a,y=1;for(let f=Math.max(t+c,i-h);f<i;f++){const p=Math.max(y,Math.trunc((i-f)*c/Math.max(1,h)));n[f]=u(f-p,f+p)}for(let f=Math.min(r+h-1,e-c-1);f>r;f--){const p=Math.max(y,Math.trunc((f-r)*c/Math.max(1,h)));n[f]=u(f-p,f+p)}for(let f=c+1;f<i-h;f++)n[f]=u(f-c,f+c);for(let f=Math.min(r+h,e-c-1);f<o-c-1;f++)n[f]=u(f-c,f+c);return n}function Us(n){return!Number.isFinite(n)||Math.abs(n)<=1e-9?1:Math.sin(n)/n}let qn=null,Vi=null;function Qo(){if(qn)return qn;const n=.625,t=1/128,e=Math.max(16,Math.round(n*2/t)+1),i=[],r=[];for(let c=0;c<e;c++){const y=-n+c*t;i.push(y),r.push(Math.abs(y)<=n?ai(y,si,.125,1):0)}const s=4,o=1/1024,a=Math.round(s/o)+1,l=new Array(a).fill(0),u=new Array(a).fill(1);let h=0;for(let c=0;c<r.length;c++)h+=r[c];h=Math.max(1e-9,h);for(let c=0;c<a;c++){const y=c*o;l[c]=y;let f=0;for(let p=0;p<i.length;p++)f+=r[p]*Math.cos(2*Math.PI*y*i[p]);u[c]=Math.max(1e-6,Math.abs(f)/h)}return qn={freqs:l,values:u},qn}function Ko(n,t){const e=Math.max(1e-6,Us(Math.PI*n*t)),i=Qo(),r=Te(Math.max(0,Math.min(i.freqs[i.freqs.length-1],n)),i.freqs,i.values);return Math.max(1e-6,e*r)}function $o(){if(Vi)return Vi;const n=new Array(ce/16*4).fill(1),t=ce*16,e=new Float32Array(t);for(let s=0;s<t;s++){const o=(s-t/2)/(16*At);e[s]=Math.abs(o)<=.625?ai(o,si,cr,1):0}const i=new gn(t);i.transform(e);const r=Math.max(1e-9,Math.abs(i._real[0]));n[0]=1;for(let s=1;s<n.length;s++){const o=Us(Math.PI*s/256),a=Math.max(1e-6,Math.hypot(i._real[s],i._imag[s])/r);n[s]=Math.max(1e-6,o*a)}return Vi=n,n}function Jo(n,t,e=Ft){if(n.length===0||t.length!==n.length)return null;const i=ce,r=i/2,s=cr,o=2*At,a=Math.max(1,Math.round(.5*At)),l=5,u=Math.max(0,Math.round(r-e*At)),h=Math.min(i-1,Math.round(r+e*At));if(h-u<32)return null;let c=new Array(i).fill(An),y=0,f=0,p=0,d=0,g=u,x=h,m=u,b=h,_=0;for(;;){const R=new Array(i).fill(0),Q=new Array(i).fill(0);c=new Array(i).fill(An),y=0,f=0,p=0,d=0;let G=-1,X=-1;for(let Y=0;Y<n.length;Y++){const pt=Math.trunc(n[Y]*At+r),mt=Math.max(m,pt-l),gt=Math.min(b-1,pt+l);for(let ht=mt;ht<=gt;ht++){const fe=(ht-r)*s,Wt=Math.max(0,1-Math.abs((n[Y]-fe)*1.75));Wt>0&&(Q[ht]+=t[Y]*Wt,R[ht]+=Wt)}}const st=Math.max(r-i/8,m+2*At),$=Math.min(r+i/8,b-2*At);for(let Y=Math.max(0,m-1);Y<=Math.min(i-1,b+1);Y++)R[Y]>0&&(c[Y]=Q[Y]/R[Y],Y<st&&(y+=c[Y],p++),Y>$&&(f+=c[Y],d++),G<0&&(G=Y),X=Y);if(G<0||X<0||p<=0||d<=0)return null;for(let Y=G-1;Y>=0;Y--)c[Y]=c[G];for(let Y=X+1;Y<i;Y++)c[Y]=c[X];const Z=new Array(i).fill(0);let lt=r;const ut=2*At;for(let Y=ut+1;Y<i-1-ut;Y++){let pt=0,mt=0;for(let gt=-ut;gt<=ut;gt++)mt+=c[Y+gt]*gt,pt+=gt*gt;Z[Y]=pt>0?mt/pt:0,Math.abs(Z[Y])>Math.abs(Z[lt]??0)&&Y>m+ut&&Y<b-ut-1&&(lt=Y)}if(Math.abs(lt-r)>2*At&&Math.abs(lt-r)<12*At)return null;let nt=0;for(let Y=Math.max(0,r-ut);Y<=Math.min(i-1,r+ut);Y++)Math.abs(Z[Y])>Math.abs(nt)&&(nt=Z[Y]);if(!Number.isFinite(nt)||Math.abs(nt)<=1e-9)return null;const Zt=Math.abs(nt*.001);g=m,x=b;let Yt=!1;for(let Y=r-o;Y>=m+a;Y--)if(Z[Y]*nt<0&&Math.abs(Z[Y])>Zt){let pt=0,mt=0,gt=0;for(let ht=Y;ht>=m;ht--)Z[ht]*nt<0&&(pt++,mt=Math.max(mt,Math.abs(Z[ht]))),gt++;if(pt>gt*.4&&mt/Math.abs(nt)>.25||pt>.9*gt&&gt>o){g=Math.min(Y,r-o),Yt=!0;break}}for(let Y=r+o;Y<b-a;Y++)if(Z[Y]*nt<0&&Math.abs(Z[Y])>Zt){let pt=0,mt=0,gt=0;for(let ht=Y;ht<b;ht++)Z[ht]*nt<0&&(pt++,mt=Math.max(mt,Math.abs(Z[ht]))),gt++;if(pt>gt*.4&&mt/Math.abs(nt)>.25||pt>.9*gt&&gt>o){x=Math.max(Y,r+o),Yt=!0;break}}if(Yt&&_<2){m=g,b=x,_++;continue}break}const P=Math.max(y/p,f/d),M=Math.min(y/p,f/d);let w=g,C=g,F=1/0,A=1/0;for(let R=g;R<=x;R++){const Q=(c[Math.max(0,R-2)]+c[Math.max(0,R-1)]+c[R]+c[Math.min(i-1,R+1)]+c[Math.min(i-1,R+2)])/5,G=Math.abs(Q-M-.1*(P-M)),X=Math.abs(Q-M-.9*(P-M));G<F&&(F=G,w=R),X<A&&(A=X,C=R)}if(w<C){const R=w;w=C,C=R}const k=Math.max(4,Math.abs(w-C)*s),v=Math.max(a,a+2*Math.trunc(k/Math.max(s,1e-6)));w+=v,C-=v;const S=Math.max(Math.abs(w-r),Math.abs(C-r),Math.max(o,Math.trunc(4/Math.max(s,1e-6)))),T=new Array(i).fill(0),L=new Array(i).fill(0),E=1.85,V=.5;for(let R=0;R<n.length;R++){const Q=Math.trunc(n[R]*At+r);let G=5;Math.abs(Q-r)>V*S&&(G=Math.abs(Q-r)>2*V*S?12:7);const X=Math.max(g,Q-G),st=Math.min(x-1,Q+G);if(st<r-E*S||X>r+E*S){for(let $=X;$<=st;$++)T[$]+=t[R],L[$]+=1;continue}for(let $=X;$<=st;$++){let Z=1;if(Math.abs($-r)<E*S){const lt=($-r)*s;if(Math.abs($-r)<S*V)Z=ai(n[R]-lt,si,s,1);else{const ut=(Math.abs($-r)/Math.max(1e-6,S)-V)/Math.max(1e-6,E-V),nt=1*(1-ut)+.01*ut;Z=ai(n[R]-lt,si,s,nt)}}!(Z>0)||!Number.isFinite(Z)||(T[$]+=t[R]*Z,L[$]+=Z)}}const U=new Array(i).fill(0);let N=-1,D=-1;for(let R=Math.max(0,g-1);R<=Math.min(i-1,x+1);R++)L[R]>0?(U[R]=T[R]/L[R],N<0&&(N=R),D=R):U[R]=An;if(N<0||D<0)return null;const B=3*At;let z=U[N],q=1;for(let R=N+1;R<r&&q<B;R++)U[R]!==An&&(z+=U[R],q++);z/=Math.max(1,q);let O=U[D],W=1;for(let R=D-1;R>r&&W<B;R--)U[R]!==An&&(O+=U[R],W++);O/=Math.max(1,W);for(let R=N-1;R>=0;R--)U[R]=z;for(let R=D+1;R<i;R++)U[R]=O;const j=Math.max(Math.trunc(r-E*S),g+2),K=Math.min(Math.trunc(r+E*S),x-3),et=Math.max(1,Math.trunc(2/Math.max(s,1e-6))),tt=qo(U,g,x,j,K,et),ot=new Array(i).fill(0);let rt=tt[Math.max(0,Math.min(i-1,g))]??tt[0]??0;for(let R=g;R<=x;R++){const Q=tt[R]??rt;ot[R]=(tt[Math.min(i-1,R+1)]??Q)-rt,rt=Q}return{esf:tt,lsfFull:ot}}function Zo(n){const t=new Array(n.length).fill(0);if(n.length===0)return t;t[0]=n[0];for(let e=1;e<n.length;e++){const i=ur(n[e]-n[e-1]);t[e]=t[e-1]+i}return t}function ur(n){if(!Number.isFinite(n))return n;let t=(n+Math.PI)%(2*Math.PI);return t<0&&(t+=2*Math.PI),t-Math.PI}function tl(n,t,e=0){if(n.length===0)return[];const i=Number.isFinite(e)?e:0,r=n.map((o,a)=>{const l=t[a]??0,u=-2*Math.PI*i*l;return ur(o-u)});return Zo(r).map((o,a)=>{const l=t[a]??0,u=-2*Math.PI*i*l;return o+u})}function el(n,t,e,i=.05,r=Number.POSITIVE_INFINITY){const s=Math.min(n.length,t.length);if(s<2)return null;let o=0;if(e)for(let d=1;d<s;d++){const g=e[d];Number.isFinite(g)&&g>o&&(o=g)}const a=e&&o>0?Math.max(1e-6,o*i):0;let l=0,u=0,h=0,c=0,y=0,f=0;for(let d=1;d<s;d++){const g=n[d],x=t[d];if(!Number.isFinite(g)||!Number.isFinite(x)||Math.abs(g)<=1e-12||g>r)continue;const m=e?e[d]:1;if(!Number.isFinite(m)||m<=a)continue;const b=e?m*m:1;l+=b,u+=b*g,h+=b*x,c+=b*g*g,y+=b*g*x,f++}if(f<4||l<=0)return null;const p=l*c-u*u;return Math.abs(p)<=1e-12?null:{slope:(l*y-u*h)/p,intercept:(h*c-u*y)/p,used:f,threshold:a}}function nl(n,t,e=Number.POSITIVE_INFINITY){const i=[],r=[],s=Math.min(n.length,t.length);for(let o=1;o<s;o++)Number.isFinite(n[o])&&Number.isFinite(t[o])&&Math.abs(n[o])>1e-12&&n[o]<=e&&(i.push(n[o]),r.push(t[o]));return i.length<2?{slope:0,intercept:Number.isFinite(t[0])?t[0]:0,used:i.length}:{...lr(i,r),used:i.length}}function il(n,t,e,i=.05,r=Number.POSITIVE_INFINITY,s=0){const o=Math.min(n.length,t.length);if(o<4)return null;let a=0;if(e)for(let b=1;b<o;b++){const _=e[b];Number.isFinite(_)&&_>a&&(a=_)}const l=e&&a>0?Math.max(1e-6,a*i):0,u=[];for(let b=1;b<o;b++){const _=n[b],P=t[b];if(!Number.isFinite(_)||!Number.isFinite(P)||Math.abs(_)<=1e-12||_>r)continue;const M=e?e[b]:1;!Number.isFinite(M)||M<=l||u.push({freq:_,phase:P,weight:e?M*M:1})}if(u.length<4)return null;const h=b=>{let _=0,P=0;for(const A of u){const k=A.phase+2*Math.PI*b*A.freq;_+=A.weight*Math.sin(k),P+=A.weight*Math.cos(k)}const M=Math.atan2(_,P);let w=0,C=0;const F=.65;for(const A of u){const k=A.phase+2*Math.PI*b*A.freq,v=Math.abs(ur(k-M)),S=v<=F?v*v:F*(2*v-F);w+=A.weight*S,C+=A.weight}return{score:C>0?w/C:Number.POSITIVE_INFINITY,intercept:M}},c=Number.isFinite(s)?s:0,y=Math.max(2,Math.min(8,Math.abs(c)>1e-6?4:2)),f=.02;let p=c,d=h(p);for(let b=c-y;b<=c+y+f*.5;b+=f){const _=h(b);_.score<d.score&&(d=_,p=b)}let g=p-f*2,x=p+f*2;for(let b=0;b<32;b++){const _=g+(x-g)/3,P=x-(x-g)/3,M=h(_).score,w=h(P).score;M<w?x=P:g=_}const m=(g+x)*.5;return d=h(m),{slope:-2*Math.PI*m,intercept:d.intercept,used:u.length,threshold:l}}function rl(n,t){if(Number.isFinite(t)&&t>0)return t;let e=0;for(const i of n)Number.isFinite(i)&&i>e&&(e=i);return e>0?e:Number.POSITIVE_INFINITY}function hr(n,t,e,i=Number.POSITIVE_INFINITY,r=0){const s=tl(n,t,r),o=rl(t,i),a=il(t,n,e,.05,o,r),l=a?null:el(t,s,e,.05,o),u=a||l?null:nl(t,s,o),h=(a==null?void 0:a.slope)??(l==null?void 0:l.slope)??(u==null?void 0:u.slope)??0,c=(a==null?void 0:a.intercept)??(l==null?void 0:l.intercept)??(u==null?void 0:u.intercept)??0,y=s.map((m,b)=>m-(h*(t[b]??0)+c)),f=Number.isFinite(y[0])?y[0]:0,p=t.map(m=>h*m+c+f),d=y.map(m=>m-f),g=c+f,x=Number.isFinite(h)?-h/(2*Math.PI):null;return{raw:[...n],unwrapped:s,linear:p,residual:d,fit:{groupDelayPx:x===null?null:x-r,absoluteGroupDelayPx:x,referenceDelayPx:r,slopeRadPerCycle:Number.isFinite(h)?h:null,interceptRad:Number.isFinite(g)?g:null,fitPointCount:(a==null?void 0:a.used)??(l==null?void 0:l.used)??(u==null?void 0:u.used)??0,fitWeightThreshold:(a==null?void 0:a.threshold)??(l==null?void 0:l.threshold)??0,fitDomain:"cycles-per-pixel",fitMaxFreqCyclesPerPixel:o}}}function dr(n,t,e){const i=[],r=[],s=[],o=[];for(let a=0;a<e.length;a++){const l=e[a];i.push(Te(l,t,n.raw)),r.push(Te(l,t,n.unwrapped)),s.push(Te(l,t,n.linear)),o.push(Te(l,t,n.residual))}return{ptfRaw:i,ptfUnwrapped:r,ptfLinear:s,ptfResidual:o}}function Ds(n,t){const e=n.map((a,l)=>({dist:a,value:t[l]})).filter(a=>Number.isFinite(a.dist)&&Number.isFinite(a.value)).sort((a,l)=>a.dist-l.dist);if(e.length===0)return{dists:[],vals:[]};const i=Math.max(1,Math.min(16,Math.floor(e.length*.1)));let r=0,s=0;for(let a=0;a<i;a++)r+=e[a].value,s+=e[e.length-1-a].value;r/=i,s/=i;const o=r<=s?e:e.map(a=>({dist:-a.dist,value:a.value})).sort((a,l)=>a.dist-l.dist);return{dists:o.map(a=>a.dist),vals:o.map(a=>a.value)}}function sl(n,t,e){const i=cn(e,"RAW green-site sampling"),r=t%2,s=n%2;return i==="RGGB"||i==="BGGR"?(r+s)%2!==0:i==="GBRG"||i==="GRBG"?(r+s)%2===0:!1}function al(n,t){return n+t&1?2:1}function ts(n,t){return(t&1)<<1|n&1}function bn(n,t,e){return n===void 0?0:typeof n=="number"?Number.isFinite(n)?n:0:Number.isFinite(n[ts(t,e)])?n[ts(t,e)]:0}function _t(n,t,e,i){return i!==void 0&&i!=="default"?al(n,t)===i:sl(n,t,e)}function oi(n){return n.length===0?0:n.reduce((t,e)=>t+e,0)/n.length}function Bs(n,t,e){const i=(e%t+t)%t,r=Math.floor(i),s=(r+1)%t,o=i-r,a=r<n.length?n[r]:0,l=s<n.length?n[s]:0;return a*(1-o)+l*o}function Os(n,t){const e=n.length;if(e===0)return 0;if(t<=0)return n[0];if(t>=e-1)return n[e-1];const i=Math.floor(t),r=Math.min(e-1,i+1),s=t-i;return n[i]*(1-s)+n[r]*s}function ol(n,t,e){const i=n.length,r=new Array(i).fill(0);for(let s=0;s<i;s++)r[s]=Os(n,s-e+t);return r}function li(n,t,e){const i=n.length;if(i===0)return{peakPos:0,peakIdx:0,peakVal:0};const r=Math.max(0,Math.floor(t-e)),s=Math.min(i-1,Math.ceil(t+e));let o=Math.max(0,Math.min(i-1,Math.round(t))),a=-1/0;for(let u=r;u<=s;u++){const h=Math.abs(n[u]);h>a&&(a=h,o=u)}Number.isFinite(a)||(a=Math.abs(n[o]??0));let l=o;if(o>0&&o<i-1){const u=n[o]>=0?1:-1,h=u*n[o-1],c=u*n[o],y=u*n[o+1],f=h-2*c+y;if(Number.isFinite(f)&&Math.abs(f)>1e-9){const p=.5*(h-y)/f;Number.isFinite(p)&&Math.abs(p)<=1&&(l=o+p)}}return{peakPos:l,peakIdx:o,peakVal:Math.abs(Os(n,l))}}function ll(n,t,e,i,r,s,o,a,l){const u=Math.floor(i.x),h=Math.floor(i.y),c=Math.floor(i.w),y=Math.floor(i.h),f=[],p=(d,g)=>{if(d<0||g<0||d>=t||g>=e)return null;const x=r+d,m=s+g;return Math.max(0,n[g*t+d]-bn(l,x,m))};for(let d=0;d<y;d++){const g=[],x=h+d;for(let m=0;m<c;m++){const b=u+m,_=r+b,P=s+x;if(_t(_,P,o,a)){g.push(p(b,x)??0);continue}const M=[],w=p(b-1,x),C=p(b+1,x),F=p(b,x-1),A=p(b,x+1);if(w!==null&&_t(_-1,P,o,a)&&M.push(w),C!==null&&_t(_+1,P,o,a)&&M.push(C),F!==null&&_t(_,P-1,o,a)&&M.push(F),A!==null&&_t(_,P+1,o,a)&&M.push(A),M.length===0){const k=[],v=p(b-1,x-1),S=p(b+1,x-1),T=p(b-1,x+1),L=p(b+1,x+1);v!==null&&_t(_-1,P-1,o,a)&&k.push(v),S!==null&&_t(_+1,P-1,o,a)&&k.push(S),T!==null&&_t(_-1,P+1,o,a)&&k.push(T),L!==null&&_t(_+1,P+1,o,a)&&k.push(L),g.push(oi(k));continue}g.push(oi(M))}f.push(g)}return f}function zs(n,t,e,i,r,s,o,a,l){const u=Math.floor(i.x),h=Math.floor(i.y),c=Math.floor(i.w),y=Math.floor(i.h),f=[];for(let p=0;p<y;p++){const d=h+p;for(let g=0;g<c;g++){const x=u+g,m=r+x,b=s+d;_t(m,b,o,a)&&f.push({x:m,y:b,value:Math.max(0,n[d*t+x]-bn(l,m,b))})}}return f}function Vs(n,t,e,i,r){const s=Math.floor(i.x),o=Math.floor(i.y),a=Math.floor(i.w),l=Math.floor(i.h),u=(r==null?void 0:r.globalX)??0,h=(r==null?void 0:r.globalY)??0,c=!!(r!=null&&r.isThreePlane)&&n.length>=t*e*3,y=r==null?void 0:r.threePlaneChannel,f=[];for(let p=0;p<l;p++){const d=o+p,g=d*t;for(let x=0;x<a;x++){const m=s+x;let b=0;if(!c)b=Math.max(0,n[g+m]-bn(r==null?void 0:r.blackLevel,u+m,h+d));else{const _=(g+m)*3;if(y!==void 0)b=n[_+y];else{const P=n[_],M=n[_+1],w=n[_+2];b=.2126*P+.7152*M+.0722*w}}f.push({x:u+m,y:h+d,value:b})}}return f}function cl(n){var s;const t=n.length,e=((s=n[0])==null?void 0:s.length)??0;let i=0,r=0;for(let o=1;o<t-1;o++)for(let a=1;a<e-1;a++)i+=Math.abs(n[o][a+1]-n[o][a-1]),r+=Math.abs(n[o+1][a]-n[o-1][a]);return{gx:i,gy:r}}function Gs(n,t,e,i,r,s,o){var y;const a=n.length,l=((y=n[0])==null?void 0:y.length)??0,u=(f,p,d,g,x,m,b)=>{const _=b?a:l,P=Math.max(0,p-3),M=Math.min(_,p+4);let w=0,C=0;for(let A=P;A<M;A++)w+=f[A],C+=A*f[A];if(w<=0)return null;const F=C/w;return b?{x:t+g*x,y:d+F*m,weight:w}:{x:d+F*m,y:e+g*x,weight:w}},h=(f,p,d,g,x,m,b)=>{const _=Math.max(3,Math.min(Math.max(3,Math.floor(d/3)),Math.max(4,Math.round(d*.12)))),P=f.map(S=>{let T=-1/0,L=-1;for(let E=0;E<S.length;E++)S[E]>T&&(T=S[E],L=E);return{peakValue:T,peakIndex:L}}),M=(p-1)*.5,w=P.map((S,T)=>({...S,index:T})).filter(S=>S.peakValue>1&&S.peakIndex>=0).sort((S,T)=>{const L=T.peakValue-S.peakValue;return Math.abs(L)>1e-6?L:Math.abs(S.index-M)-Math.abs(T.index-M)});if(w.length===0)return[];const C=w[0],F=new Array(p).fill(null),A=u(f[C.index],C.peakIndex,g,C.index,x,m,b);if(!A)return[];F[C.index]=A;const k=(S,T)=>{const L=f[S],E=P[S];if(!(E.peakValue>1)||E.peakIndex<0)return null;const V=Math.max(0,Math.floor(T-_)),U=Math.min(L.length,Math.ceil(T+_+1));let N=-1/0,D=-1;for(let O=V;O<U;O++)L[O]>N&&(N=L[O],D=O);if(D<0||!(N>1))return null;const B=Math.max(1e-6,E.peakValue),z=Math.abs(D-T)<=_,q=N>=B*.25;return!z||!q?null:u(L,D,g,S,x,m,b)};let v=A?b?(A.y-g)/m:(A.x-g)/m:C.peakIndex;for(let S=C.index+1;S<p;S++){const T=k(S,v);T&&(F[S]=T,v=b?(T.y-g)/m:(T.x-g)/m)}v=A?b?(A.y-g)/m:(A.x-g)/m:C.peakIndex;for(let S=C.index-1;S>=0;S--){const T=k(S,v);T&&(F[S]=T,v=b?(T.y-g)/m:(T.x-g)/m)}return F.filter(S=>!!S)};if(s){const f=n.map(p=>p.map((d,g)=>g===0?0:Math.abs(d-p[g-1])));return h(f,a,l,t,r,i,!1)}const c=[];for(let f=0;f<l;f++){const p=new Array(a).fill(0);for(let d=1;d<a;d++)p[d]=Math.abs(n[d][f]-n[d-1][f]);c.push(p)}return h(c,l,a,e,i,r,!0)}function be(n){if(n.length<2)return null;let t=0,e=0,i=0;for(const c of n)t+=c.weight,e+=c.x*c.weight,i+=c.y*c.weight;if(t<=0)return null;e/=t,i/=t;let r=0,s=0,o=0;for(const c of n){const y=c.x-e,f=c.y-i;r+=c.weight*y*y,s+=c.weight*f*f,o+=c.weight*y*f}r/=t,s/=t,o/=t;const a=.5*Math.atan2(2*o,r-s);let l=Math.cos(a),u=Math.sin(a);const h=Math.hypot(l,u);return!Number.isFinite(h)||h<=1e-9?null:(l/=h,u/=h,(l<0||Math.abs(l)<=1e-9&&u<0)&&(l=-l,u=-u),{pointX:e,pointY:i,dirX:l,dirY:u,orientation:Math.abs(l)>=Math.abs(u)?1:2})}function ul(n,t){if(n.length!==4||t.length!==4||n.some(i=>i.length!==4))return null;const e=n.map((i,r)=>[...i,t[r]]);for(let i=0;i<4;i++){let r=i,s=Math.abs(e[i][i]);for(let a=i+1;a<4;a++){const l=Math.abs(e[a][i]);l>s&&(s=l,r=a)}if(!(s>1e-12))return null;if(r!==i){const a=e[i];e[i]=e[r],e[r]=a}const o=e[i][i];for(let a=i;a<=4;a++)e[i][a]/=o;for(let a=0;a<4;a++){if(a===i)continue;const l=e[a][i];if(!(Math.abs(l)<=1e-12))for(let u=i;u<=4;u++)e[a][u]-=l*e[i][u]}}return[e[0][4],e[1][4],e[2][4],e[3][4]]}function hl(n){if(n.length<4)return 0;const t=[...n].sort((y,f)=>y.x-f.x),e=t[0].x,r=t[t.length-1].x-e;if(!(r>1e-6))return 0;const s=16,o=[];for(let y=0;y<s;y++){const f=Math.max(0,Math.floor((y-1.5)*t.length/s)),p=Math.min(t.length-1,Math.floor((y+2.5)*t.length/s));if(p<f)continue;let d=0,g=0,x=0;for(let m=f;m<=p;m++)d+=t[m].x,g+=t[m].y,x++;x>0&&o.push({x:d/x,y:g/x})}if(o.length<4)return 0;const a=[.05952381,0,-.03571429,-.04761905,-.03571429,0,.05952381],l=new Array(o.length).fill(0),u=3;for(let y=0;y<o.length;y++){let f=0;for(let p=-u;p<=u;p++){const d=y+p;d<0||d>=o.length||(f+=a[p+u]*o[d].y)}l[y]=f}let h=0,c=1/0;for(let y=0;y<l.length-1;y++){const f=l[y],p=l[y+1];if(f===0){const b=Math.abs(o[y].x);b<c&&(c=b,h=o[y].x);continue}if(f*p>=0)continue;const d=p-f,g=Math.abs(d)>1e-12?-f/d:.5,x=o[y].x+(o[y+1].x-o[y].x)*g,m=Math.abs(x);m<c&&(c=m,h=x)}return!Number.isFinite(h)||h<e+.3*r||h>e+.7*r?0:h}function dl(n){if(n.length<8)return null;const t=[...n].filter(p=>Number.isFinite(p.x)&&Number.isFinite(p.y)&&Number.isFinite(p.weight)).sort((p,d)=>p.x-d.x);if(t.length<8)return null;const i=[hl(t),0,.5*(t[Math.floor((t.length-1)*.5)].x+t[Math.ceil((t.length-1)*.5)].x)];let r=null;for(const p of i){if(!Number.isFinite(p))continue;let d=0,g=0;for(const x of t)x.x<=p?d++:g++;if(d>=4&&g>=4){r=p;break}}if(r===null)return null;const s=Array.from({length:4},()=>new Array(4).fill(0)),o=new Array(4).fill(0);for(const p of t){const d=p.x,g=p.y,x=Math.max(1e-6,p.weight),m=d<=r?[d*d,d,1,0]:[2*r*d-r*r,d,1,(d-r)*(d-r)];for(let b=0;b<4;b++){o[b]+=x*m[b]*g;for(let _=0;_<4;_++)s[b][_]+=x*m[b]*m[_]}}const a=ul(s,o);if(!a)return null;const[l,u,h,c]=a,y=u+2*(l-c)*r,f=h+(c-l)*r*r;return[l,u,h,c,y,f].every(p=>Number.isFinite(p))?{splitX:r,left:[l,u,h],right:[c,y,f]}:null}function fl(n,t,e){const[i,r,s]=e;if(Math.abs(i)<=1e-12){const b=1+r*r;return b>1e-12?[(n-r*(s-t))/b]:[n]}const o=2*i*i,a=3*i*r,l=1+2*i*s-2*i*t+r*r,u=r*s-t*r-n;if(Math.abs(o)<=1e-12)return[n];const h=a/o,c=l/o,y=u/o,f=(h*h-3*c)/9,p=(2*h*h*h-9*h*c+27*y)/54,d=p*p-f*f*f;if(d<0&&f>0){const b=Math.acos(Math.max(-1,Math.min(1,p/Math.sqrt(f*f*f)))),_=-2*Math.sqrt(f);return[_*Math.cos(b/3)-h/3,_*Math.cos((b+2*Math.PI)/3)-h/3,_*Math.cos((b-2*Math.PI)/3)-h/3]}const g=Math.sqrt(Math.max(0,d)),x=-Math.sign(p||1)*Math.cbrt(Math.abs(p)+g),m=Math.abs(x)<=1e-12?0:f/x;return[x+m-h/3]}function Xs(n,t){if(n.length<8)return null;const e=-t.dirY,i=t.dirX,r=n.map(o=>({x:(o.x-t.pointX)*t.dirX+(o.y-t.pointY)*t.dirY,y:(o.x-t.pointX)*e+(o.y-t.pointY)*i,weight:o.weight})),s=dl(r);return s?{...t,normalX:e,normalY:i,splitX:s.splitX,left:s.left,right:s.right}:null}function pl(n,t){const e=n.x-t.pointX,i=n.y-t.pointY,r=e*t.dirX+i*t.dirY,s=e*t.normalX+i*t.normalY,o=r<t.splitX?t.left:t.right,a=fl(r,s,o);let l=s,u=Number.POSITIVE_INFINITY;for(const h of a){if(!Number.isFinite(h))continue;const c=o[0]*h*h+o[1]*h+o[2],y=r-h,f=s-c,p=Math.hypot(y,f);Number.isFinite(p)&&p<u&&(u=p,l=(f>=0?1:-1)*p)}return Number.isFinite(u)?l:s}function Ys(n,t,e,i,r,s,o,a=Ft){if(!t||t.length<8||n.length===0)return null;const l=t.filter(S=>Number.isFinite(S.x)&&Number.isFinite(S.y)).map(S=>({x:S.x,y:S.y,weight:1}));if(l.length<8)return null;const u=be(l);if(!u)return null;const h=e.p2.x-e.p1.x,c=e.p2.y-e.p1.y,y=Math.hypot(h,c);if(!Number.isFinite(y)||y<=1e-6)return null;let f=u.dirX,p=u.dirY;f*h+p*c<0&&(f=-f,p=-p);const d={...u,dirX:f,dirY:p},g=Xs(l,d),x=h/y,m=c/y,b=-m,_=x,P=(e.p1.x+e.p2.x)*.5,M=(e.p1.y+e.p2.y)*.5,w=-d.dirY,C=d.dirX,F=Math.abs(x)>=Math.abs(m)?1:2,A=[],k=[];for(const S of n){const T=S.x-P,L=S.y-M,E=T*x+L*m;if(Math.abs(E)>i)continue;const V=T*b+L*_;if(Math.abs(V)>r)continue;const U=g?pl(S,g):(S.x-d.pointX)*w+(S.y-d.pointY)*C;Number.isFinite(U)&&(A.push(U),k.push(S.value))}if(A.length<8)return null;const v=o?s!=null&&s.forceLegacyModel?wn(A,k,F,a,r*2):_n(A,k,F,a):gi(A,k,Math.max(2,r*2),s==null?void 0:s.manualBinSize,F,s==null?void 0:s.preferAutoPerEdgeBin);return v?(v.quadraticProjectionUsed=!!g,v):null}function es(n){if(n.length<2)return null;const t=n.filter(e=>Number.isFinite(e.x)&&Number.isFinite(e.y)).map(e=>({x:e.x,y:e.y,weight:1}));return t.length<2?null:Tn(t,be(t))}function Tn(n,t){if(!t||n.length<2)return null;let e=1/0,i=-1/0;for(const s of n){const o=(s.x-t.pointX)*t.dirX+(s.y-t.pointY)*t.dirY;e=Math.min(e,o),i=Math.max(i,o)}if(!Number.isFinite(e)||!Number.isFinite(i))return null;const r=Math.max(.5,(i-e)*.03);return{p1:{x:t.pointX+t.dirX*(e-r),y:t.pointY+t.dirY*(e-r)},p2:{x:t.pointX+t.dirX*(i+r),y:t.pointY+t.dirY*(i+r)}}}function ml(n,t,e,i,r,s,o,a){return[Tn(n,t),Tn(e,i),Tn(r,s),Tn(o,a)]}function gl(n,t,e){if(!n||n.length<8)return;const i=n.filter(m=>Number.isFinite(m.x)&&Number.isFinite(m.y)).map(m=>({x:m.x,y:m.y,weight:1}));if(i.length<8)return;const r=be(i);if(!r)return;const s=t.p2.x-t.p1.x,o=t.p2.y-t.p1.y,a=Math.hypot(s,o);if(!Number.isFinite(a)||a<=1e-6)return;let l=r.dirX,u=r.dirY;l*s+u*o<0&&(l=-l,u=-u);const h=Xs(i,{...r,dirX:l,dirY:u});if(!h)return;const c=n.map(m=>(m.x-h.pointX)*h.dirX+(m.y-h.pointY)*h.dirY).filter(m=>Number.isFinite(m)),y=(t.p1.x-h.pointX)*h.dirX+(t.p1.y-h.pointY)*h.dirY,f=(t.p2.x-h.pointX)*h.dirX+(t.p2.y-h.pointY)*h.dirY;if(Number.isFinite(y)&&c.push(y),Number.isFinite(f)&&c.push(f),c.length<2)return;const p=Math.min(...c),d=Math.max(...c);if(!Number.isFinite(p)||!Number.isFinite(d)||d-p<=1e-6)return;const g=Math.max(21,e),x=[];for(let m=0;m<g;m++){const b=g===1?.5:m/(g-1),_=p+(d-p)*b,P=_<h.splitX?h.left:h.right,M=P[0]*_*_+P[1]*_+P[2];x.push({x:h.pointX+_*h.dirX+M*h.normalX,y:h.pointY+_*h.dirY+M*h.normalY})}return x}function yl(n,t,e,i){const r=Es(n,t);if(!r||r.profilePositionsPx.length<2||r.lsf.length<8)return null;let s=0,o=Number.POSITIVE_INFINITY;for(let a=0;a<r.profilePositionsPx.length;a++){const l=Math.abs(r.profilePositionsPx[a]);l<o&&(o=l,s=a)}return{esf:r.esf,lsfFull:r.lsf,binSize:r.diagnostics.profileStepPx,orientation:e,zeroIndex:s,shortSidePx:i,fallbackUsed:!1,continuousMtfV2:!0,continuousProjectedDistsPx:[...n],continuousProjectedValues:[...t],continuousMtfV2Result:r,mtfEngine:"continuous-v2"}}function Ws(n,t,e,i,r,s){if(n==="continuous-v2"){const a=yl(t,e,i,r);if(a)return a}const o=s();return o&&(o.mtfEngine="legacy"),o}function gi(n,t,e,i,r,s=!1,o=!1,a=!1){return Ws(a?"legacy":"continuous-v2",n,t,r,e,()=>xl(n,t,e,i,r,s,o))}function xl(n,t,e,i,r,s=!1,o=!1){if(n.length===0||t.length!==n.length)return null;const a=Ds(n,t),l=a.dists,u=a.vals;if(l.length===0)return null;const h=()=>{const P=e/2;let M=0;for(const C of n)Math.abs(C)<=P&&M++;if(M<=0)return .125;const w=40*P/M;return Math.max(.01,Math.min(.125,w))},c=(P,M,w,C,F)=>{if(!(C>0)||!(F>0)||!(w>M))return!1;const A=Math.floor((w-M)/C);if(A<2)return!1;const k=Math.max(M,-F),v=Math.min(w,F);if(!(v>k))return!1;const S=Math.max(0,Math.floor((k-M)/C)),T=Math.min(A,Math.ceil((v-M)/C));if(T<=S)return!1;const L=new Array(T-S).fill(0),E=M+S*C,V=M+T*C;for(let U=0;U<P.length;U++){const N=P[U];if(N<E)continue;if(N>=V)break;const D=Math.floor((N-M)/C);D>=S&&D<T&&L[D-S]++}return L.every(U=>U>0)},y=()=>{const P=l[0],M=l[l.length-1],w=Math.max(0,e*.25),C=.125,F=.5,A=.001,k=Math.round((F-C)/A);for(let v=0;v<=k;v++){const S=Number((C+v*A).toFixed(3));if(c(l,P,M,S,w))return S}return F};let f=.125;i&&i>0?f=Math.max(.01,Math.min(.5,i)):s?f=y():f=h();const p=l[0],d=l[l.length-1],g=Math.floor((d-p)/f);if(g<2)return null;const m=(()=>{const P=new Array(g).fill(0),M=new Array(g).fill(0);for(let C=0;C<l.length;C++){const F=(l[C]-p)/f;if(Number.isFinite(F))if(o){const A=Math.floor(F),k=F-A,v=1-k,S=k;A>=0&&A<g&&(P[A]+=u[C]*v,M[A]+=v);const T=A+1;T>=0&&T<g&&(P[T]+=u[C]*S,M[T]+=S)}else{const A=Math.floor(F);A>=0&&A<g&&(P[A]+=u[C],M[A]++)}}let w=u[0];for(let C=0;C<g;C++)M[C]>0?(P[C]/=M[C],w=P[C]):P[C]=w;return P})(),b=jo(m,f),_=Math.max(0,Math.min(g-1,-p/f-.5));return{esf:m,lsfFull:b,binSize:f,orientation:r,zeroIndex:_,shortSidePx:e,fallbackUsed:!0}}function _n(n,t,e,i=Ft){if(n.length===0||t.length!==n.length)return null;const r=Ds(n,t),s=r.dists.map((l,u)=>({dist:l,value:r.vals[u]})).filter(l=>Math.abs(l.dist)<i);if(s.length<8)return null;const o=s.map(l=>l.dist),a=s.map(l=>l.value);return Ws("continuous-v2",o,a,e,i*2,()=>wn(o,a,e,i))}function wn(n,t,e,i=Ft,r=i*2){if(n.length===0||t.length!==n.length)return null;const s=[],o=[];for(let a=0;a<n.length;a++){const l=n[a],u=t[a];!Number.isFinite(l)||!Number.isFinite(u)||Math.abs(l)>=i||(s.push(l),o.push(u))}return s.length<8?null:gi(s,o,Math.max(2,r),void 0,e,!0,!0,!0)}function xe(n,t,e,i,r,s,o,a,l=Ft){if(s<=0||o<=0)return null;const h=!(!!(a!=null&&a.isThreePlane)&&n.length>=t*e*3)&&((a==null?void 0:a.greenOnly)??!1),c=h?cn(a==null?void 0:a.bayerPattern,"constrained RAW SFR sampling"):null,y=r.p2.x-r.p1.x,f=r.p2.y-r.p1.y,p=Math.hypot(y,f);if(!Number.isFinite(p)||p<=1e-6)return null;const d=y/p,g=f/p,x=-g,m=d,b=(r.p1.x+r.p2.x)*.5,_=(r.p1.y+r.p2.y)*.5,P=Math.abs(d)>=Math.abs(g)?1:2,M=h?zs(n,t,e,i,0,0,c,a==null?void 0:a.greenPhase,a==null?void 0:a.blackLevel):Vs(n,t,e,i,{...a,globalX:0,globalY:0});if(M.length===0)return null;if(!(a!=null&&a.disableQuadraticProjection)){const F=Ys(M,a==null?void 0:a.quadraticFitPoints,r,s,o,a,!0,l);if(F)return F}const w=[],C=[];for(const F of M){const A=F.x-b,k=F.y-_,v=A*d+k*g;if(Math.abs(v)>s)continue;const S=A*x+k*m;Math.abs(S)>o||(w.push(S),C.push(F.value))}return w.length<8?null:a!=null&&a.forceLegacyModel?wn(w,C,P,l):_n(w,C,P,l)}function bl(n,t,e=0){const i=[...n.lsfFull];if(i.length<3)return!1;const r=Math.max(n.binSize,1e-6),s=Number.isFinite(n.zeroIndex)?n.zeroIndex:i.length/2,o=Math.max(1,Math.round((n.shortSidePx??0)*.5/r));let{peakPos:a,peakIdx:l,peakVal:u}=li(i,s,o);const h=u*.2;let c=0,y=i.length-1;for(let d=l;d>=0;d--)if(i[d]<h){c=d;break}for(let d=l;d<i.length;d++)if(i[d]<h){y=d;break}const f=y-c;if(t&&f>0){const d=f*4,g=[],x=[];if(e>0){const m=Math.max(0,l-d-e),b=Math.max(0,l-d);for(let M=m;M<b;M++)g.push(M),x.push(i[M]);const _=Math.min(i.length,l+d),P=Math.min(i.length,l+d+e);for(let M=_;M<P;M++)g.push(M),x.push(i[M])}else{for(let m=0;m<Math.max(0,l-d);m++)g.push(m),x.push(i[m]);for(let m=Math.min(i.length,l+d);m<i.length;m++)g.push(m),x.push(i[m])}if(g.length>2){const{slope:m,intercept:b}=lr(g,x);for(let _=0;_<i.length;_++)i[_]=i[_]-(m*_+b);({peakPos:a}=li(i,s,o))}}return Math.abs(a-s)*r<=Math.max(1e-6,(n.shortSidePx??0)/6)}function _l(n){const t=n.length;if(t<3)return!1;let e=0,i=-1/0;for(let o=0;o<t;o++){const a=Math.abs(n[o]);a>i&&(i=a,e=o)}const r=t/3,s=2*t/3;return e>=r&&e<=s}function xt(n,t,e){const i=Math.max(0,Math.floor(n.x)),r=Math.max(0,Math.floor(n.y)),s=Math.min(t,Math.ceil(n.x+n.w)),o=Math.min(e,Math.ceil(n.y+n.h)),a=s-i,l=o-r;return a<2||l<2?null:{x:i,y:r,w:a,h:l}}function fr(n,t,e,i){const r=[],s=n.x,o=n.y,a=n.x+n.w,l=n.y+n.h,u=n.x+n.w*.5,h=n.y+n.h*.5,c=[{x:s,y:o},{x:a,y:o},{x:a,y:l},{x:s,y:l},{x:u,y:o},{x:a,y:h},{x:u,y:l},{x:s,y:h},{x:u,y:h}];for(const y of c){const f=we(y,t);Number.isFinite(f.x)&&Number.isFinite(f.y)&&r.push(f)}return r.length===0?null:xt(Nt(r,2),e,i)}function Nt(n,t=0){let e=1/0,i=1/0,r=-1/0,s=-1/0;for(const o of n)e=Math.min(e,o.x),i=Math.min(i,o.y),r=Math.max(r,o.x),s=Math.max(s,o.y);return{x:e-t,y:i-t,w:r-e+t*2,h:s-i+t*2}}function ns(n,t){let e=Math.atan2(t,n)*180/Math.PI;return e<0&&(e+=180),e}function he(n,t){const e=n.p2.x-n.p1.x,i=n.p2.y-n.p1.y,r=Math.hypot(e,i);if(!Number.isFinite(r)||r<=1e-6)return null;const s=-i/r,o=e/r;return[{x:n.p1.x+s*t,y:n.p1.y+o*t},{x:n.p2.x+s*t,y:n.p2.y+o*t},{x:n.p2.x-s*t,y:n.p2.y-o*t},{x:n.p1.x-s*t,y:n.p1.y-o*t}]}function wl(n,t,e,i,r,s,o,a){if(s<=0||o<=0)return null;const u=!(!!(a!=null&&a.isThreePlane)&&n.length>=t*e*3)&&((a==null?void 0:a.greenOnly)??!1),h=u?cn(a==null?void 0:a.bayerPattern,"constrained RAW edge sampling"):null,c=r.p2.x-r.p1.x,y=r.p2.y-r.p1.y,f=Math.hypot(c,y);if(!Number.isFinite(f)||f<=1e-6)return null;const p=c/f,d=y/f,g=-d,x=p,m=(r.p1.x+r.p2.x)*.5,b=(r.p1.y+r.p2.y)*.5,_=Math.abs(p)>=Math.abs(d)?1:2,P=u?zs(n,t,e,i,0,0,h,a==null?void 0:a.greenPhase,a==null?void 0:a.blackLevel):Vs(n,t,e,i,{...a,globalX:0,globalY:0});if(P.length===0)return null;if(!(a!=null&&a.disableQuadraticProjection)){const C=Ys(P,a==null?void 0:a.quadraticFitPoints,r,s,o,a,!1);if(C)return C}const M=[],w=[];for(const C of P){const F=C.x-m,A=C.y-b,k=F*p+A*d;if(Math.abs(k)>s)continue;const v=F*g+A*x;Math.abs(v)>o||(M.push(v),w.push(C.value))}return M.length<8?null:gi(M,w,o*2,a==null?void 0:a.manualBinSize,_,a==null?void 0:a.preferAutoPerEdgeBin)}function Ml(n,t,e){const i=[...n],r=new Array(n.length).fill(0),s=[0,0,0];let o=-1,a=1,l=-1;for(let c=1;c<n.length-1;c++){let y=0;if(n[c]>1e-4){y=Math.atan2(e[c]*o,t[c]*o);let f=0;for(let p=-5;p<=5;p++)Math.abs(y+p*2*Math.PI-s[1])<Math.abs(y+f*2*Math.PI-s[1])&&(f=p);y+=f*2*Math.PI}c>3&&Math.abs(y-s[0])>Math.PI/2&&l<c-1&&n[c]<.5&&(a*=-1,l=c),i[c]*=a,o*=-1,s[0]=s[1],s[1]=y,s[2]=y}const u=[-.086,.343,.486,.343,-.086];for(let c=0;c<n.length-3;c++){let y=0;for(let f=-2;f<=2;f++)y+=i[Math.abs(c+f)]*u[f+2];r[c]=y}for(let c=0;c<n.length-3;c++)i[c]=r[c];const h=7;for(let c=0;c<3;c++){r.fill(0);for(let f=0;f<n.length-h;f++)if(f<h)r[f]=i[f];else{const p=Math.min(5,Math.floor((f-5)/3)),d=Ho[p];let g=0;for(let x=-h;x<=h;x++)g+=i[f+x]*d[x+h];r[f]=g}for(let f=n.length-h-2;f<n.length;f++)r[f]=i[f];const y=Math.abs(r[0])>1e-9?r[0]:1;for(let f=0;f<n.length;f++)i[f]=r[f]/y}for(let c=0;c<n.length;c++)i[c]=Math.abs(i[c]);return i}function Sl(n,t){const e=[[0,0,0],[0,0,0],[0,0,0]],i=[0,0,0];for(let o=0;o<n.length;o++){const a=n[o],l=-t+o,u=[1,a,a*a];for(let h=0;h<3;h++){i[h]+=u[h]*l;for(let c=0;c<3;c++)e[h][c]+=u[h]*u[c]}}const r=aa(e);if(!r)return null;const s=oa(r,i);return[s[0],s[1],s[2]]}function vl(n,t){let e=0,i=1,r=0,s=!1,o=0;const a=Math.min(n.length,ce/16*2);for(let l=0;l<a&&!s;l++){const u=n[l];if(i>.5&&u<=.5){const h=-(u-i)*ce;Math.abs(h)>1e-9&&(r=-((.5-i-h*e)/h),o=l,s=!0)}i=u,e=l/ce}if(!s)return null;if(o>=5&&o<a-10){const l=Math.min(Math.max(2,o-9),9),u=Sl(n.slice(o-l,o+l+1),l);if(u){const c=(u[0]+.5*u[1]+.25*u[2]+o)/ce;if(o>9)r=c;else{const f=(o-5)/ce/8;r=(1-f)*r+f*c}}}return r*At*t}function Pl(n,t,e){for(let i=1;i<t.length;i++){if(t[i-1]<e||t[i]>=e)continue;const r=t[i-1]-t[i];if(Math.abs(r)<=1e-12)return n[i];const s=(t[i-1]-e)/r;return n[i-1]+s*(n[i]-n[i-1])}return null}function Cl(n,t){const i=new gn(4096),r=new Float32Array(4096);let s=0;for(let d=1;d<n.lsf.length;d++)Math.abs(n.lsf[d])>Math.abs(n.lsf[s])&&(s=d);for(let d=0;d<4096;d++)r[d]=Bs(n.lsf,4096,d+s);i.transform(r);const o=Math.max(1e-12,Math.hypot(i._real[0],i._imag[0])),a=Math.max(1e-9,n.diagnostics.profileStepPx),l=Math.max(...t),u=Math.min(4096/2,Math.ceil(l*4096*a)+2),h=[],c=[],y=[];for(let d=0;d<=u;d++)h.push(d/(4096*a)),c.push(d===0?0:Math.atan2(i._imag[d],i._real[d])),y.push(Math.hypot(i._real[d],i._imag[d])/o);const f=hr(c,h,y,l,0);return{...dr(f,h,t),fit:f.fit}}function kl(n,t){const e=[];for(const u of n){const h=u.continuousMtfV2Result;if(h){e.push(h);continue}if(!u.continuousProjectedDistsPx||!u.continuousProjectedValues)continue;const c=Es(u.continuousProjectedDistsPx,u.continuousProjectedValues);c&&e.push(c)}if(e.length===0)return null;const i=new Array(501),r=new Array(501).fill(0);for(let u=0;u<=500;u++){const h=u/250;i[u]=h;for(const c of e)r[u]+=Te(h,c.frequencies,c.mtf);r[u]/=e.length}const s=1,o=e.length===1?e[0].mtf50:Pl(i,r,.5),a=e[0],l=Cl(a,i);return{esf:a.esf,lsf:[],lsfCropped:a.lsf,mtf:r,ptf:l.ptfResidual,ptfRaw:l.ptfRaw,ptfUnwrapped:l.ptfUnwrapped,ptfLinear:l.ptfLinear,ptfResidual:l.ptfResidual,ptfPhaseFit:l.fit,freqs:i.map(u=>u*s),mtf50:o===null?null:o*s,calcRadius:e.reduce((u,h)=>u+h.diagnostics.halfWidthPx,0)/e.length}}function Fl(n,t){if(n.length===0)return null;const e=ce,i=ce/16*4,r=new gn(e),s=1,o=$o(),a=new Float32Array(501);for(let v=0;v<=500;v++)a[v]=v/500*s*2;const l=new Array(i).fill(0).map((v,S)=>S/e*s*At),u=new Float32Array(i).fill(0),h=new Float32Array(i).fill(0);let c=0,y=[],f=[],p=[],d=[],g=[],x=[],m=[],b=null,_=0;for(const v of n){const S=v.mtfmapperOrderedDists&&v.mtfmapperOrderedVals&&v.mtfmapperOrderedDists.length===v.mtfmapperOrderedVals.length?Jo(v.mtfmapperOrderedDists,v.mtfmapperOrderedVals,v.mtfmapperEffectiveMaxDot??Ft):null,T=(S==null?void 0:S.lsfFull)??v.lsfFull,L=(S==null?void 0:S.esf)??v.esf;if(T.length<e)continue;const E=new Float32Array(e);for(let N=0;N<e;N++)E[N]=T[N]??0;r.transform(E);const V=Math.max(1e-9,Math.abs(r._real[0])),U=new Array(i).fill(0);for(let N=1;N<i;N++)U[N]=Math.atan2(r._imag[N],r._real[N]);for(let N=0;N<i;N++)u[N]+=r._real[N]/V,h[N]+=r._imag[N]/V;if(c++,_+=v.shortSidePx*.5,y.length===0){y=[...T],f=[...L];const N=new Array(i).fill(0);N[0]=1;for(let O=1;O<i;O++)N[O]=Math.hypot(r._real[O]/V,r._imag[O]/V);const D=l.map(O=>O),B=(Number.isFinite(v.zeroIndex)?v.zeroIndex:0)*(v.binSize??cr),z=hr(U,D,N,Number.POSITIVE_INFINITY,B),q=dr(z,l,a);d=q.ptfRaw,g=q.ptfUnwrapped,x=q.ptfLinear,m=q.ptfResidual,p=q.ptfResidual,b=z.fit}}if(c===0)return null;const P=new Float32Array(i),M=new Float32Array(i),w=new Array(i).fill(0);w[0]=1;for(let v=0;v<i;v++)P[v]=u[v]/c,M[v]=h[v]/c,v>0&&(w[v]=Math.hypot(P[v],M[v]));const C=Ml(w,P,M),F=new Array(i).fill(0);for(let v=0;v<i;v++)F[v]=C[v]/o[v];const A=Array.from(a,v=>Te(v,l,F)),k=vl(F,s);return{esf:f,lsf:[],lsfCropped:y,mtf:A,ptf:p,ptfRaw:d,ptfUnwrapped:g,ptfLinear:x,ptfResidual:m,ptfPhaseFit:b,freqs:Array.from(a),mtf50:k,calcRadius:_/c}}function Al(n,t,e,i=!1,r=0,s=!1){if(n.length===0)return null;if(n.every(w=>w.mtfEngine==="continuous-v2"||w.continuousMtfV2))return kl(n);if(n.every(w=>w.mtfmapperLike))return Fl(n);const o=4096,a=new gn(o),l=1,u=new Float32Array(501);for(let w=0;w<=500;w++)u[w]=w/500*l*2;const h=new Float32Array(501).fill(0);let c=0,y=[],f=[],p=0,d=[],g=[],x=[],m=[],b=[],_=null;for(const w of n){let C=[...w.lsfFull];const F=w.binSize,A=Number.isFinite(w.zeroIndex)?w.zeroIndex:C.length/2,k=Math.max(1,Math.round((w.shortSidePx??0)*.5/Math.max(F,1e-6)));let{peakPos:v,peakIdx:S,peakVal:T}=li(C,A,k);const L=T*.2;let E=0,V=C.length-1;for(let R=S;R>=0;R--)if(C[R]<L){E=R;break}for(let R=S;R<C.length;R++)if(C[R]<L){V=R;break}const U=V-E;let N=!1;if(i&&U>0){const R=U*4,Q=[],G=[];if(r>0){const X=Math.max(0,S-R-r),st=Math.max(0,S-R);for(let lt=X;lt<st;lt++)Q.push(lt),G.push(C[lt]);const $=Math.min(C.length,S+R),Z=Math.min(C.length,S+R+r);for(let lt=$;lt<Z;lt++)Q.push(lt),G.push(C[lt])}else{for(let X=0;X<Math.max(0,S-R);X++)Q.push(X),G.push(C[X]);for(let X=Math.min(C.length,S+R);X<C.length;X++)Q.push(X),G.push(C[X])}if(Q.length>2){const{slope:X,intercept:st}=lr(Q,G);for(let $=0;$<C.length;$++)C[$]=C[$]-(X*$+st);({peakPos:v,peakIdx:S,peakVal:T}=li(C,A,k)),N=!0}}let D=0,B=0;if(t>0)B=t,D=Math.round(t/F);else{const R=T*.2;let Q=0,G=C.length-1;for(let Z=S;Z>=0;Z--)if(C[Z]<R){Q=Z;break}for(let Z=S;Z<C.length;Z++)if(C[Z]<R){G=Z;break}const st=(G-Q)*F;let $=Math.max(2,st*8);B=$,D=Math.round($/F)}p+=B;const z=Math.max(0,Math.floor(A-D)),q=Math.min(C.length,Math.ceil(A+D)),O=C.slice(z,q);if(O.length<8)continue;const W=new Float32Array(o).fill(0),j=new Array(O.length).fill(0);for(let R=0;R<O.length;R++){let Q=1;s&&(Q=.5*(1-Math.cos(2*Math.PI*R/(O.length-1)))),j[R]=O[R]*Q}const K=Math.max(0,Math.min(O.length-1,v-z));for(let R=0;R<o;R++)W[R]=Bs(j,o,R+K);a.transform(W);const et=[],tt=[],ot=[];for(let R=0;R<=o/2;R++){const Q=a._real[R],G=a._imag[R],X=Math.sqrt(Q*Q+G*G);et.push(X),tt.push(R/(o*F)*l),ot.push(Math.atan2(G,Q))}const rt=et[0];if(rt>0){for(let R=0;R<=500;R++){const Q=u[R],X=Ko(Q,F);h[R]+=Te(Q,tt,et)/rt/X}if(c++,y.length===0){y=ol(O,K,(O.length-1)/2),f=N?Tl(C):w.esf;const R=et.map(st=>st/rt),Q=tt.map(st=>st),G=hr(ot,Q,R,Number.POSITIVE_INFINITY,0),X=dr(G,tt,u);g=X.ptfRaw,x=X.ptfUnwrapped,m=X.ptfLinear,b=X.ptfResidual,d=X.ptfResidual,_=G.fit}}}if(c===0)return null;const P=Array.from(h).map(w=>w/c);let M=null;for(let w=0;w<P.length-1;w++)if(P[w]>=.5&&P[w+1]<.5){M=u[w]+(.5-P[w])*(u[w+1]-u[w])/(P[w+1]-P[w]);break}return{esf:f,lsf:[],lsfCropped:y,mtf:P,ptf:d,ptfRaw:g,ptfUnwrapped:x,ptfLinear:m,ptfResidual:b,ptfPhaseFit:_,freqs:Array.from(u),mtf50:M,calcRadius:p/c}}function Tl(n){const t=new Array(n.length).fill(0);let e=0;for(let i=0;i<n.length;i++)e+=n[i],t[i]=e;return t}function Te(n,t,e){if(n<=t[0])return e[0];if(n>=t[t.length-1])return e[e.length-1];let i=0;for(;n>t[i+1];)i++;const r=(n-t[i])/(t[i+1]-t[i]);return e[i]+r*(e[i+1]-e[i])}function pr(n){return{...Qi,...n,gradientPercentiles:n!=null&&n.gradientPercentiles&&n.gradientPercentiles.length>0?n.gradientPercentiles:Qi.gradientPercentiles}}function Il(n){return!n||n.length===0?void 0:[Number.isFinite(n[0])?n[0]:0,Number.isFinite(n[1])?n[1]:Number.isFinite(n[0])?n[0]:0,Number.isFinite(n[2])?n[2]:Number.isFinite(n[0])?n[0]:0,Number.isFinite(n[3])?n[3]:Number.isFinite(n[0])?n[0]:0]}function Nl(n,t){const e=n.width,i=n.height,r=n.data,s=cn(n.bayerPattern,"RAW SFR detection"),o=Il(n.blackLevels),a=new Float32Array(e*i),l=(_,P)=>_<0||P<0||_>=e||P>=i?null:Math.max(0,r[P*e+_]-bn(o,_,P));let u=1/0,h=-1/0;for(let _=0;_<i;_++){const P=_*e;for(let M=0;M<e;M++){const w=P+M;let C=0;if(_t(M,_,s,t))C=l(M,_)??0;else{const F=[],A=l(M-1,_),k=l(M+1,_),v=l(M,_-1),S=l(M,_+1);if(A!==null&&_t(M-1,_,s,t)&&F.push(A),k!==null&&_t(M+1,_,s,t)&&F.push(k),v!==null&&_t(M,_-1,s,t)&&F.push(v),S!==null&&_t(M,_+1,s,t)&&F.push(S),F.length>0)C=oi(F);else{const T=[],L=l(M-1,_-1),E=l(M+1,_-1),V=l(M-1,_+1),U=l(M+1,_+1);L!==null&&_t(M-1,_-1,s,t)&&T.push(L),E!==null&&_t(M+1,_-1,s,t)&&T.push(E),V!==null&&_t(M-1,_+1,s,t)&&T.push(V),U!==null&&_t(M+1,_+1,s,t)&&T.push(U),C=oi(T)}}a[w]=C,C<u&&(u=C),C>h&&(h=C)}}if(!Number.isFinite(u)||!Number.isFinite(h)||h<=u+1e-9)return new Uint8Array(e*i);const c=1024,y=new Uint32Array(c),f=h-u;for(let _=0;_<a.length;_++){const P=Math.max(0,Math.min(1,(a[_]-u)/f)),M=Math.min(c-1,Math.max(0,Math.floor(P*(c-1))));y[M]++}const p=a.length,d=_=>{const P=p*_;let M=0;for(let w=0;w<c;w++)if(M+=y[w],M>=P)return u+w/Math.max(1,c-1)*f;return h},g=d(.01),x=d(.99),m=Math.max(1e-9,x-g),b=new Uint8Array(e*i);for(let _=0;_<a.length;_++){const P=Math.max(0,Math.min(1,(a[_]-g)/m));b[_]=Math.round(P*255)}return b}function Rl(n,t,e){const i=new Float32Array(n.length),r=new Float32Array(n.length),s=new Float32Array(n.length);for(let o=1;o<e-1;o++)for(let a=1;a<t-1;a++){const l=o*t+a,u=n[(o-1)*t+(a-1)],h=n[(o-1)*t+a],c=n[(o-1)*t+(a+1)],y=n[o*t+(a-1)],f=n[o*t+(a+1)],p=n[(o+1)*t+(a-1)],d=n[(o+1)*t+a],g=n[(o+1)*t+(a+1)],x=-u-2*y-p+(c+2*f+g),m=-u-2*h-c+(p+2*d+g);i[l]=x,r[l]=m,s[l]=Math.hypot(x,m)}return{gx:i,gy:r,magnitude:s}}function Ll(n,t){let e=0,i=0;for(let l=0;l<n.length;l++){const u=n[l];!Number.isFinite(u)||u<=1e-6||(e=Math.max(e,u),i++)}if(i===0||e<=1e-6)return[];const r=1024,s=new Uint32Array(r);for(let l=0;l<n.length;l++){const u=n[l];if(!Number.isFinite(u)||u<=1e-6)continue;const h=Math.max(0,Math.min(1,u/e)),c=Math.min(r-1,Math.floor(h*(r-1)));s[c]++}const o=t&&t.length>0?t:Qi.gradientPercentiles,a=[];for(const l of o){const u=i*l;let h=0;for(let c=0;c<r;c++)if(h+=s[c],h>=u){a.push(c/Math.max(1,r-1)*e);break}}return Array.from(new Set(a.filter(l=>l>0))).sort((l,u)=>u-l)}function El(n,t){const e=new Uint8Array(n.length);for(let i=0;i<n.length;i++)e[i]=n[i]>=t?1:0;return e}const Ul=256*256;function Dl(n,t,e){if(n.length>=Ul){const s=vo.compute(n,t,e);if(s)return{gray:s.blurredGray,gradient:{gx:s.gx,gy:s.gy,magnitude:s.magnitude},backend:"webgl"}}const i=Hl(n,t,e),r=Rl(i,t,e);return{gray:i,gradient:r,backend:"cpu"}}function Bl(n,t,e,i){let r=n;for(let s=0;s<i;s++){const o=new Uint8Array(n.length);for(let a=0;a<e;a++)for(let l=0;l<t;l++){let u=0;for(let h=-1;h<=1&&!u;h++){const c=a+h;if(!(c<0||c>=e))for(let y=-1;y<=1;y++){const f=l+y;if(!(f<0||f>=t)&&r[c*t+f]){u=1;break}}}o[a*t+l]=u}r=o}return r}function Ol(n,t,e){const i=new Int32Array(n.length),r=[];let s=1;for(let o=0;o<n.length;o++){if(!n[o]||i[o]!==0)continue;const a=[o];i[o]=s;let l=0,u=t,h=e,c=0,y=0,f=0,p=!1;for(;l<a.length;){const d=a[l++],g=d%t,x=Math.floor(d/t);f++,u=Math.min(u,g),h=Math.min(h,x),c=Math.max(c,g),y=Math.max(y,x),(g===0||x===0||g===t-1||x===e-1)&&(p=!0);for(let m=-1;m<=1;m++)for(let b=-1;b<=1;b++){if(b===0&&m===0)continue;const _=g+b,P=x+m;if(_<0||P<0||_>=t||P>=e)continue;const M=P*t+_;!n[M]||i[M]!==0||(i[M]=s,a.push(M))}}r.push({label:s,x:u,y:h,w:c-u+1,h:y-h+1,area:f,touchesBorder:p}),s++}return{labels:i,components:r}}function js(n,t){const e=Math.hypot(n,t);if(!Number.isFinite(e)||e<=1e-9)return null;let i=n/e,r=t/e;return(i<0||Math.abs(i)<=1e-9&&r<0)&&(i=-i,r=-r),{x:i,y:r}}function ve(n,t){if(n.length===0)return 0;const e=[...n].sort((o,a)=>o.value-a.value),i=e.reduce((o,a)=>o+Math.max(0,a.weight),0);if(i<=0)return e[Math.floor((e.length-1)*t)].value;const r=Math.max(0,Math.min(1,t))*i;let s=0;for(const o of e)if(s+=Math.max(0,o.weight),s>=r)return o.value;return e[e.length-1].value}function is(n){const t=n.filter(i=>Number.isFinite(i)).sort((i,r)=>i-r);if(t.length===0)return 0;const e=i=>{if(i.length===1)return i[0];if(i.length===2)return(i[0]+i[1])*.5;const r=Math.ceil(i.length*.5);let s=0,o=1/0;for(let a=0;a+r-1<i.length;a++){const l=i[a+r-1]-i[a];l<o&&(o=l,s=a)}return e(i.slice(s,s+r))};return e(t)}function zl(n,t,e,i,r,s,o){const a=[];for(let l=r.y;l<r.y+r.h;l++)for(let u=r.x;u<r.x+r.w;u++){const h=l*s+u;if(n[h]!==t||!e[h])continue;const c=i.magnitude[h];!Number.isFinite(c)||c<=1e-6||a.push({x:u,y:l,weight:c,gx:i.gx[h],gy:i.gy[h]})}return a}function Vl(n){let t=0,e=0,i=0,r=0,s=0;for(const l of n){t+=l.weight,e+=l.x*l.weight,i+=l.y*l.weight;const u=Math.hypot(l.gx,l.gy);if(!Number.isFinite(u)||u<=1e-6)continue;const h=-l.gy/u,c=l.gx/u;r+=l.weight*(h*h-c*c),s+=l.weight*(2*h*c)}if(t<=0)return null;e/=t,i/=t;const o=.5*Math.atan2(s,r),a=js(Math.cos(o),Math.sin(o));return a?{centerX:e,centerY:i,dirX:a.x,dirY:a.y,orthoX:-a.y,orthoY:a.x}:null}function Qn(n,t){let e=0,i=0;const r=-t.dirY,s=t.dirX;for(const o of n){const a=(o.x-t.pointX)*r+(o.y-t.pointY)*s;i+=o.weight*a*a,e+=o.weight}return e<=0?1/0:Math.sqrt(i/e)}function rs(n,t,e,i,r){const s=Math.max(0,Math.min(t-1,i)),o=Math.max(0,Math.min(e-1,r)),a=Math.floor(s),l=Math.floor(o),u=Math.min(t-1,a+1),h=Math.min(e-1,l+1),c=s-a,y=o-l,f=n[l*t+a],p=n[l*t+u],d=n[h*t+a],g=n[h*t+u],x=f+(p-f)*c,m=d+(g-d)*c;return x+(m-x)*y}function Gl(n,t,e,i,r,s,o,a){const l=rs(n,t,e,i-s*a,r-o*a);return rs(n,t,e,i+s*a,r+o*a)-l}function Kn(n,t,e,i,r){const s=Math.max(1e-6,e-t);if(n.length===0||!Number.isFinite(s))return{points:[],coverageRatio:0,centerCoverageRatio:0};const o=Math.max(1.5,Math.min(4,s/18)),a=Math.max(1,Math.ceil(s/o)),l=new Map;for(const p of n){const d=i(p);if(!Number.isFinite(d)||d<t||d>e)continue;const g=Math.max(0,Math.min(a-1,Math.floor((d-t)/o))),x=p.weight/(1+Math.abs(r(p))),m=l.get(g);(!m||x>m.score)&&l.set(g,{point:p,score:x})}const u=Array.from(l.values()).sort((p,d)=>i(p.point)-i(d.point)).map(p=>p.point),h=Math.max(0,Math.floor(a*.3)),c=Math.max(h+1,Math.ceil(a*.7));let y=0;for(let p=h;p<c;p++)l.has(p)&&y++;const f=Math.max(1,c-h);return{points:u,coverageRatio:u.length/a,centerCoverageRatio:y/f}}function $n(n,t){const e=n.dirX*t.dirY-n.dirY*t.dirX;if(!Number.isFinite(e)||Math.abs(e)<=1e-6)return null;const i=t.pointX-n.pointX,r=t.pointY-n.pointY,s=(i*t.dirY-r*t.dirX)/e;return{x:n.pointX+n.dirX*s,y:n.pointY+n.dirY*s}}function Xl(n){if(n.length<3)return 0;let t=0;for(let e=0;e<n.length;e++){const i=n[e],r=n[(e+1)%n.length];t+=i.x*r.y-r.x*i.y}return t*.5}function Yl(n,t,e,i,r,s,o,a,l){const u=zl(i,r.label,s,o,r,t),h=u.map(I=>({x:I.x,y:I.y}));if(u.length<l.minEdgePoints)return{candidate:null,failureStage:"min_edge_points",pointsCount:u.length,strongEdgePoints:h};const c=Vl(u);if(!c)return{candidate:null,failureStage:"dominant_axes",pointsCount:u.length,strongEdgePoints:h};const y=u.map(I=>{const it=I.x-c.centerX,at=I.y-c.centerY;return{...I,u:it*c.dirX+at*c.dirY,v:it*c.orthoX+at*c.orthoY}}),f={x:c.centerX,y:c.centerY},p=ve(y.map(I=>({value:I.u,weight:I.weight})),l.extentQuantileLow),d=ve(y.map(I=>({value:I.u,weight:I.weight})),l.extentQuantileHigh),g=ve(y.map(I=>({value:I.v,weight:I.weight})),l.extentQuantileLow),x=ve(y.map(I=>({value:I.v,weight:I.weight})),l.extentQuantileHigh),m=Math.max(1e-6,Math.max(Math.abs(p),Math.abs(d))),b=Math.max(1e-6,Math.max(Math.abs(g),Math.abs(x))),_=72,P=360/_,M=Array.from({length:_},()=>[]),w=I=>{let it=I%360;return it<0&&(it+=360),it},C=(I,it)=>{const at=Math.abs(w(I)-w(it));return Math.min(at,360-at)};y.forEach(I=>{const it=I.u/m,at=I.v/b,vt=w(Math.atan2(at,it)*180/Math.PI),Et=Math.hypot(it,at),kt=Math.max(0,Math.min(_-1,Math.floor(vt/P)));M[kt].push({point:I,angleDeg:vt,normRadius:Et})});const F=M.map(I=>I.length>0?is(I.map(it=>it.normRadius)):-1/0),A=(I,it)=>{let at=-1,vt=-1/0;for(let ft=0;ft<M.length;ft++){if(M[ft].length===0)continue;const se=(ft+.5)*P;if(C(se,I)>45||it.some(ia=>C(se,ia)<45))continue;const Se=F[ft];Se>vt&&(vt=Se,at=ft)}let Et=at>=0?(at+.5)*P:I,kt=at>=0?M[at]:y.map(ft=>{const re=ft.u/m,se=ft.v/b;return{point:ft,angleDeg:w(Math.atan2(se,re)*180/Math.PI),normRadius:Math.hypot(re,se)}}).filter(ft=>C(ft.angleDeg,I)<=45&&!it.some(re=>C(ft.angleDeg,re)<45));if(kt.length===0&&(kt=y.map(ft=>{const re=ft.u/m,se=ft.v/b;return{point:ft,angleDeg:w(Math.atan2(se,re)*180/Math.PI),normRadius:Math.hypot(re,se)}}).filter(ft=>C(ft.angleDeg,I)<=45),Et=I),kt.length===0)return{x:y[0].x,y:y[0].y,u:y[0].u,v:y[0].v,angleDeg:I};const Yn=at>=0?F[at]:is(kt.map(ft=>ft.normRadius));let ge=0,Cn=0,Sr=0,vr=0,Pr=0;for(const ft of kt){const re=C(ft.angleDeg,I)/45,se=Math.abs(ft.normRadius-Yn),Se=Math.max(1e-6,ft.point.weight)/(1+re*2+se*6);ge+=Se,Cn+=ft.point.x*Se,Sr+=ft.point.y*Se,vr+=ft.point.u*Se,Pr+=ft.point.v*Se}return ge>0?{x:Cn/ge,y:Sr/ge,u:vr/ge,v:Pr/ge,angleDeg:Et}:{x:kt[0].point.x,y:kt[0].point.y,u:kt[0].point.u,v:kt[0].point.v,angleDeg:kt[0].angleDeg}},k=A(225,[]),v=A(315,[k.angleDeg]),S=A(45,[k.angleDeg,v.angleDeg]),T=A(135,[k.angleDeg,v.angleDeg,S.angleDeg]),L=[{x:k.x,y:k.y},{x:v.x,y:v.y},{x:S.x,y:S.y},{x:T.x,y:T.y}],E=d-p,V=x-g,U=Math.min(E,V),N=Math.max(E,V);if(!Number.isFinite(U)||U<l.minSpanPx||N/Math.max(1,U)>l.maxAspectRatio)return{candidate:null,failureStage:"span_aspect",pointsCount:u.length,minSpan:U,maxSpan:N,axisCentroid:f,axisExtremePoints:L,strongEdgePoints:h};const D=Math.max(l.bandMinPx,Math.min(l.bandMaxPx,U*l.bandScale)),B=Math.max(1,Math.min(3,D*.55)),z=Math.max(a,0),q=void 0,O=void 0,W=I=>I.map(it=>({x:it.x,y:it.y,weight:it.weight})),j=I=>I.map(it=>({x:it.x,y:it.y})),K=(I,it,at)=>I.filter(vt=>{if(!Number.isFinite(vt.weight)||vt.weight<z)return!1;const Et=Gl(n,t,e,vt.x,vt.y,it,at,B);return Number.isFinite(Et)&&Et>=l.minPointContrast}),et=p,tt=d,ot=g,rt=x,R=l.minCoverageRatio,Q=l.minCenterCoverageRatio,G=[],X=[],st=[],$=[],Z=[],lt=(I,it,at,vt,Et,kt)=>(at-I)*(kt-it)-(vt-it)*(Et-I),ut=I=>I>1e-6?1:I<-1e-6?-1:0,nt=[{u:(k.u+v.u)*.5,v:(k.v+v.v)*.5},{u:(v.u+S.u)*.5,v:(v.v+S.v)*.5},{u:(S.u+T.u)*.5,v:(S.v+T.v)*.5},{u:(T.u+k.u)*.5,v:(T.v+k.v)*.5}],Zt=(I,it)=>{const at=ut(lt(k.u,k.v,S.u,S.v,I,it)),vt=ut(lt(v.u,v.v,T.u,T.v,I,it));return`${at},${vt}`},Yt=new Map;nt.forEach((I,it)=>{Yt.set(Zt(I.u,I.v),it)});for(const I of y){if(!Number.isFinite(I.u)||!Number.isFinite(I.v)){Z.push(I);continue}let at=Yt.get(Zt(I.u,I.v))??-1;if(at<0){let vt=1/0;for(let Et=0;Et<nt.length;Et++){const kt=nt[Et],Yn=(I.u-kt.u)/m,ge=(I.v-kt.v)/b,Cn=Yn*Yn+ge*ge;Cn<vt&&(vt=Cn,at=Et)}}at===0?G.push(I):at===1?X.push(I):at===2?st.push(I):at===3?$.push(I):Z.push(I)}const Y=[...G,...st],pt=[...X,...$],mt={dir:Y.length,ortho:pt.length,unassigned:y.length-Y.length-pt.length},gt=G.length>=l.minSidePoints?ve(G.map(I=>({value:I.v,weight:I.weight})),.5):g,ht=st.length>=l.minSidePoints?ve(st.map(I=>({value:I.v,weight:I.weight})),.5):x,fe=$.length>=l.minSidePoints?ve($.map(I=>({value:I.u,weight:I.weight})),.5):p,Wt=X.length>=l.minSidePoints?ve(X.map(I=>({value:I.u,weight:I.weight})),.5):d,Ot=[{x:(k.x+v.x)*.5,y:(k.y+v.y)*.5},{x:(v.x+S.x)*.5,y:(v.y+S.y)*.5},{x:(S.x+T.x)*.5,y:(S.y+T.y)*.5},{x:(T.x+k.x)*.5,y:(T.y+k.y)*.5}],Re=G.filter(I=>Math.abs(I.v-gt)<=D&&I.u>=et&&I.u<=tt),Le=st.filter(I=>Math.abs(I.v-ht)<=D&&I.u>=et&&I.u<=tt),Mn=$.filter(I=>Math.abs(I.u-fe)<=D&&I.v>=ot&&I.v<=rt),un=X.filter(I=>Math.abs(I.u-Wt)<=D&&I.v>=ot&&I.v<=rt),Rn=[Re.length,un.length,Le.length,Mn.length],Sn=[j(Re),j(un),j(Le),j(Mn)],Ee=K(Re,-c.orthoX,-c.orthoY),Ue=K(Le,c.orthoX,c.orthoY),pe=K(Mn,-c.dirX,-c.dirY),Me=K(un,c.dirX,c.dirY),vn=[Ee.length,Me.length,Ue.length,pe.length],De=[j(Ee),j(Me),j(Ue),j(pe)],Be=Kn(Ee,et,tt,I=>I.u,I=>I.v-gt),Oe=Kn(Me,ot,rt,I=>I.v,I=>I.u-Wt),ze=Kn(Ue,et,tt,I=>I.u,I=>I.v-ht),Ve=Kn(pe,ot,rt,I=>I.v,I=>I.u-fe),Ln=(I,it)=>I.slice().sort((at,vt)=>it(at)-it(vt)),En=Ln(Ee,I=>I.u),Un=Ln(Me,I=>I.v),Dn=Ln(Ue,I=>I.u),Bn=Ln(pe,I=>I.v),yi=[En.length,Un.length,Dn.length,Bn.length],zt=[Be.coverageRatio,Oe.coverageRatio,ze.coverageRatio,Ve.coverageRatio];Be.centerCoverageRatio,Oe.centerCoverageRatio,ze.centerCoverageRatio,Ve.centerCoverageRatio;const Rt=[j(En),j(Un),j(Dn),j(Bn)],jt={axisPointCounts:mt,sideBandPointCounts:Rn,sideContrastPointCounts:vn,gradientThreshold:a,pointAxisMinDot:q,pointAxisMargin:O,bandWidth:D,minPointContrast:l.minPointContrast,minCoverageRatio:R,minCenterCoverageRatio:Q,axisCentroid:f,axisExtremePoints:L,axisSideCenters:Ot,strongEdgePoints:h,axisDirPoints:j(Y),axisOrthoPoints:j(pt),axisUnassignedPoints:j(Z),sideBandPoints:Sn,sideContrastPoints:De};if(En.length<l.minSidePoints||Dn.length<l.minSidePoints||Bn.length<l.minSidePoints||Un.length<l.minSidePoints)return{candidate:null,failureStage:"min_side_points",pointsCount:u.length,minSpan:U,maxSpan:N,sidePointCounts:yi,sideCoverageRatios:zt,...jt,sideFitPoints:Rt};if(Be.coverageRatio<R||Oe.coverageRatio<R||ze.coverageRatio<R||Ve.coverageRatio<R||Be.centerCoverageRatio<Q||Oe.centerCoverageRatio<Q||ze.centerCoverageRatio<Q||Ve.centerCoverageRatio<Q)return{candidate:null,failureStage:"side_coverage",pointsCount:u.length,minSpan:U,maxSpan:N,sidePointCounts:yi,sideCoverageRatios:zt,...jt,sideFitPoints:Rt};const xi=En,bi=Dn,_i=Bn,wi=Un,qt=yi,Ge=be(W(xi)),Xe=be(W(bi)),Ye=be(W(_i)),We=be(W(wi));if(!Ge||!Xe||!Ye||!We)return{candidate:null,failureStage:"fit_lines",pointsCount:u.length,minSpan:U,maxSpan:N,sidePointCounts:qt,sideCoverageRatios:zt,...jt,sideFitPoints:Rt};const me=ml(W(xi),Ge,W(wi),We,W(bi),Xe,W(_i),Ye),On=l.minAxisDot,zn=(I,it,at)=>Math.abs(I.dirX*it+I.dirY*at),Lt=[zn(Ge,c.dirX,c.dirY),zn(We,c.orthoX,c.orthoY),zn(Xe,c.dirX,c.dirY),zn(Ye,c.orthoX,c.orthoY)];if(Lt[0]<On||Lt[1]<On||Lt[2]<On||Lt[3]<On)return{candidate:null,failureStage:"axis_alignment",pointsCount:u.length,minSpan:U,maxSpan:N,sidePointCounts:qt,sideCoverageRatios:zt,axisDots:Lt,...jt,sideFitPoints:Rt,sideFitLines:me};const Qt=Math.max(l.residualLimitFloor,D*l.residualLimitScale),dt=[Qn(W(xi),Ge),Qn(W(bi),Xe),Qn(W(_i),Ye),Qn(W(wi),We)],Mi=[dt[0],dt[3],dt[1],dt[2]],Si=Math.max(...dt),te=$n(Ge,Ye),ee=$n(Ge,We),ne=$n(Xe,We),ie=$n(Xe,Ye);if(!te||!ee||!ne||!ie)return{candidate:null,failureStage:"corners",pointsCount:u.length,minSpan:U,maxSpan:N,sidePointCounts:qt,sideCoverageRatios:zt,axisDots:Lt,sideResiduals:[dt[0],dt[3],dt[1],dt[2]],residualLimit:Qt,...jt,sideFitPoints:Rt,sideFitLines:me};const Pn=[te,ee,ne,ie],Kt=Math.abs(Xl(Pn));if(!Number.isFinite(Kt)||Kt<l.minQuadArea)return{candidate:null,failureStage:"quad_area",pointsCount:u.length,minSpan:U,maxSpan:N,sidePointCounts:qt,sideCoverageRatios:zt,axisDots:Lt,sideResiduals:[dt[0],dt[3],dt[1],dt[2]],residualLimit:Qt,...jt,sideFitPoints:Rt,sideFitLines:me,quadArea:Kt};const vi=Math.hypot(ee.x-te.x,ee.y-te.y),Pi=Math.hypot(ne.x-ee.x,ne.y-ee.y),Ci=Math.hypot(ne.x-ie.x,ne.y-ie.y),ki=Math.hypot(ie.x-te.x,ie.y-te.y),je=[vi,Pi,Ci,ki],Fi=Math.min(vi,Pi,Ci,ki),Zs=Math.max(vi,Pi,Ci,ki);if(!Number.isFinite(Fi)||Fi<l.minSideLength||Zs/Math.max(1,Fi)>l.maxAspectRatio)return{candidate:null,failureStage:"side_length",pointsCount:u.length,minSpan:U,maxSpan:N,sidePointCounts:qt,sideCoverageRatios:zt,axisDots:Lt,sideResiduals:[dt[0],dt[3],dt[1],dt[2]],residualLimit:Qt,...jt,sideFitPoints:Rt,sideFitLines:me,quadArea:Kt,sideLengths:je};const He=js(ee.x-te.x+(ne.x-ie.x),ee.y-te.y+(ne.y-ie.y));if(!He)return{candidate:null,failureStage:"corners",pointsCount:u.length,minSpan:U,maxSpan:N,sidePointCounts:qt,sideCoverageRatios:zt,axisDots:Lt,sideResiduals:[dt[0],dt[3],dt[1],dt[2]],residualLimit:Qt,...jt,sideFitPoints:Rt,sideFitLines:me,quadArea:Kt,sideLengths:je};const _r={x:-He.y,y:He.x},Ai=(te.x+ee.x+ne.x+ie.x)*.25,Ti=(te.y+ee.y+ne.y+ie.y)*.25,Vn=Pn.map(I=>{const it=I.x-Ai,at=I.y-Ti;return{u:it*He.x+at*He.y,v:it*_r.x+at*_r.y}}),Gn=(Math.max(...Vn.map(I=>I.u))-Math.min(...Vn.map(I=>I.u)))*.5,Xn=(Math.max(...Vn.map(I=>I.v))-Math.min(...Vn.map(I=>I.v)))*.5;if(!Number.isFinite(Gn)||!Number.isFinite(Xn)||Math.min(Gn,Xn)<6)return{candidate:null,failureStage:"box_size",pointsCount:u.length,minSpan:U,maxSpan:N,sidePointCounts:qt,sideCoverageRatios:zt,axisDots:Lt,sideResiduals:[dt[0],dt[3],dt[1],dt[2]],residualLimit:Qt,sideFitPoints:Rt,quadArea:Kt,sideLengths:je};const H=$l(n,t,e,Pn,Ge,We,Xe,Ye,Ai,Ti,Gn,Xn,l.innerPurityStdScale,l.outerMeanSpreadLimit);if(!Number.isFinite(Si)||Si>Qt)return{candidate:null,failureStage:"residual",pointsCount:u.length,minSpan:U,maxSpan:N,sidePointCounts:qt,sideCoverageRatios:zt,axisDots:Lt,sideResiduals:Mi,residualLimit:Qt,...jt,sideFitPoints:Rt,sideFitLines:me,quadArea:Kt,sideLengths:je,outerContrast:H.contrast,outerUniformityOk:H.ok,outerMeanSpread:H.meanSpread,outerMeanSpreadLimit:H.meanSpreadLimit,outerAvgStd:H.avgStd,outerAvgStdLimit:H.avgStdLimit,outerSideMeans:H.outerSideMeans,outerSideStds:H.outerSideStds,outerSideStdLimit:H.outerSideStdLimit,outerSideQuads:H.outerSideQuads,innerSideUniformityOk:H.innerSideOk,innerSideStds:H.innerSideStds,innerSideStdLimit:H.innerSideStdLimit,innerSideQuads:H.innerSideQuads};const wr=l.filterBlockPurity&&(!H.ok||!H.innerSideOk);if(wr||H.contrast<l.minOuterContrast)return{candidate:null,failureStage:wr?H.ok?"inner_roi_uniformity":"outer_uniformity":"outer_contrast",pointsCount:u.length,minSpan:U,maxSpan:N,sidePointCounts:qt,sideCoverageRatios:zt,axisDots:Lt,sideResiduals:Mi,residualLimit:Qt,...jt,sideFitPoints:Rt,sideFitLines:me,quadArea:Kt,sideLengths:je,outerContrast:H.contrast,outerUniformityOk:H.ok,outerMeanSpread:H.meanSpread,outerMeanSpreadLimit:H.meanSpreadLimit,outerAvgStd:H.avgStd,outerAvgStdLimit:H.avgStdLimit,outerSideMeans:H.outerSideMeans,outerSideStds:H.outerSideStds,outerSideStdLimit:H.outerSideStdLimit,outerSideQuads:H.outerSideQuads,innerSideUniformityOk:H.innerSideOk,innerSideStds:H.innerSideStds,innerSideStdLimit:H.innerSideStdLimit,innerSideQuads:H.innerSideQuads};const Mr=xt(Nt(Pn,1),t,e);if(!Mr)return{candidate:null,failureStage:"bbox",pointsCount:u.length,minSpan:U,maxSpan:N,sidePointCounts:qt,sideCoverageRatios:zt,axisDots:Lt,sideResiduals:[dt[0],dt[3],dt[1],dt[2]],residualLimit:Qt,...jt,sideFitPoints:Rt,sideFitLines:me,quadArea:Kt,sideLengths:je,outerContrast:H.contrast,outerUniformityOk:H.ok,outerMeanSpread:H.meanSpread,outerMeanSpreadLimit:H.meanSpreadLimit,outerAvgStd:H.avgStd,outerAvgStdLimit:H.avgStdLimit,outerSideMeans:H.outerSideMeans,outerSideStds:H.outerSideStds,outerSideStdLimit:H.outerSideStdLimit,outerSideQuads:H.outerSideQuads,innerSideUniformityOk:H.innerSideOk,innerSideStds:H.innerSideStds,innerSideStdLimit:H.innerSideStdLimit,innerSideQuads:H.innerSideQuads};const ta=1/(1+Si/Math.max(1,Qt)),ea=l.filterBlockPurity?H.score:1,na=H.contrast*ea*ta*Math.sqrt(Kt);return{candidate:{centerX:Ai,centerY:Ti,dirX:He.x,dirY:He.y,halfWidth:Gn,halfHeight:Xn,score:na,bbox:Mr,corners:Pn,sideFitPoints:Rt,outerSideMeans:H.outerSideMeans,outerSideQuads:H.outerSideQuads},failureStage:null,pointsCount:u.length,minSpan:U,maxSpan:N,sidePointCounts:qt,sideCoverageRatios:zt,axisDots:Lt,sideResiduals:Mi,residualLimit:Qt,...jt,sideFitPoints:Rt,sideFitLines:me,quadArea:Kt,sideLengths:je,outerContrast:H.contrast,outerUniformityOk:H.ok,outerMeanSpread:H.meanSpread,outerMeanSpreadLimit:H.meanSpreadLimit,outerAvgStd:H.avgStd,outerAvgStdLimit:H.avgStdLimit,outerSideMeans:H.outerSideMeans,outerSideStds:H.outerSideStds,outerSideStdLimit:H.outerSideStdLimit,outerSideQuads:H.outerSideQuads,innerSideUniformityOk:H.innerSideOk,innerSideStds:H.innerSideStds,innerSideStdLimit:H.innerSideStdLimit,innerSideQuads:H.innerSideQuads}}function Wl(n,t,e,i,r,s,o,a,l){return Yl(n,t,e,i,r,s,o,a,l).candidate}function jl(n,t,e,i,r,s){const o=pr(r),a=Math.max(i*8,i+128),l=ql(n,t,e,o.downsampleMaxSide);s==null||s("Detecting candidates: downsampling...",.02),s==null||s("Detecting candidates: edge stage...",.06);const u=Dl(l.gray,l.width,l.height),h=u.gray,c=u.gradient;s==null||s(`Detecting candidates: gradient (${u.backend==="webgl"?"WebGL1":"CPU"})...`,.1);const y=Ll(c.magnitude,o.gradientPercentiles),f=l.width*l.height,p=Math.max(o.minComponentAreaPx,Math.round(f*o.minComponentAreaRatio)),d=Math.max(p+1,Math.round(f*o.maxComponentAreaRatio)),g=[],x=Math.max(1,y.reduce((k,v,S)=>k+(S<=1,2),0));let m=0;for(let k=0;k<y.length;k++){const v=y[k],S=El(c.magnitude,v),T=k<=1?[3,2]:[2,1];for(const L of T){const E=m/x;s==null||s(`Detecting candidates: threshold ${k+1}/${y.length}, dilate ${L}`,.12+.78*E);const V=Bl(S,l.width,l.height,L),{labels:U,components:N}=Ol(V,l.width,l.height);for(const D of N){if(D.touchesBorder||D.area<p||D.area>d)continue;const B=Wl(h,l.width,l.height,U,D,S,c,v,o);if(!B)continue;const z=1/l.scale,q=B.corners.map(O=>({x:O.x*z,y:O.y*z}));g.push({centerX:B.centerX*z,centerY:B.centerY*z,dirX:B.dirX,dirY:B.dirY,halfWidth:B.halfWidth*z,halfHeight:B.halfHeight*z,score:B.score,bbox:{x:B.bbox.x*z,y:B.bbox.y*z,w:B.bbox.w*z,h:B.bbox.h*z},corners:q,sideFitPoints:B.sideFitPoints?[B.sideFitPoints[0].map(O=>({x:O.x*z,y:O.y*z})),B.sideFitPoints[1].map(O=>({x:O.x*z,y:O.y*z})),B.sideFitPoints[2].map(O=>({x:O.x*z,y:O.y*z})),B.sideFitPoints[3].map(O=>({x:O.x*z,y:O.y*z}))]:void 0,outerSideMeans:B.outerSideMeans,outerSideQuads:B.outerSideQuads?[B.outerSideQuads[0].map(O=>({x:O.x*z,y:O.y*z})),B.outerSideQuads[1].map(O=>({x:O.x*z,y:O.y*z})),B.outerSideQuads[2].map(O=>({x:O.x*z,y:O.y*z})),B.outerSideQuads[3].map(O=>({x:O.x*z,y:O.y*z}))]:void 0})}g.length>a&&(g.sort((D,B)=>B.score-D.score),g.length=a),m++}}console.log(`[SFR Auto Detect] Candidate pool before dedupe: ${g.length}`),s==null||s(`Detecting candidates: deduplicating (0/${Math.max(1,Math.min(g.length,Math.max(i*4,i+32)))})...`,.94),g.sort((k,v)=>v.score-k.score);const b=Math.max(i*4,i+32),_=g.length>b?g.slice(0,b):g,P=[];if(_.length<=256){console.log(`[SFR Auto Detect] Using simple dedupe for ${_.length} candidates`);for(let k=0;k<_.length;k++){const v=_[k];console.log(`[SFR Auto Detect] Simple dedupe candidate ${k+1}/${_.length}`,v.bbox);const S=_.length<=0?1:k/_.length;if(s==null||s(`Detecting candidates: deduplicating (${k}/${_.length})...`,.94+.05*Math.min(1,S)),!P.some(L=>{const E=Math.hypot(v.centerX-L.centerX,v.centerY-L.centerY),V=Math.max(Math.hypot(v.bbox.w,v.bbox.h),Math.hypot(L.bbox.w,L.bbox.h));return ss(v.bbox,L.bbox)>.28||E<V*.18})&&(P.push(v),P.length>=i))break}return s==null||s("Detecting candidates: deduplicating...",1),P}const M=Math.max(32,Math.round(Math.sqrt(Math.max(1,t*e)/4096))),w=new Map,C=new Set,F=k=>Math.floor(k/M),A=(k,v)=>{if(!Number.isFinite(k.bbox.x)||!Number.isFinite(k.bbox.y)||!Number.isFinite(k.bbox.w)||!Number.isFinite(k.bbox.h)||k.bbox.w<=0||k.bbox.h<=0||k.bbox.w>t*4||k.bbox.h>e*4)return;const S=F(k.bbox.x),T=F(k.bbox.x+k.bbox.w),L=F(k.bbox.y),E=F(k.bbox.y+k.bbox.h);for(let V=L;V<=E;V++)for(let U=S;U<=T;U++){const N=`${U},${V}`,D=w.get(N);D?D.push(v):w.set(N,[v])}};for(let k=0;k<_.length;k++){const v=_[k];if(k===0||k%200===0){const U=_.length<=0?1:k/_.length;s==null||s(`Detecting candidates: deduplicating (${k}/${_.length})...`,.94+.05*Math.min(1,U))}C.clear();const S=F(v.bbox.x),T=F(v.bbox.x+v.bbox.w),L=F(v.bbox.y),E=F(v.bbox.y+v.bbox.h);let V=!1;for(let U=L-1;U<=E+1&&!V;U++)for(let N=S-1;N<=T+1&&!V;N++){const D=w.get(`${N},${U}`);if(D)for(const B of D){if(C.has(B))continue;C.add(B);const z=P[B];if(!z)continue;const q=Math.hypot(v.centerX-z.centerX,v.centerY-z.centerY),O=Math.max(Math.hypot(v.bbox.w,v.bbox.h),Math.hypot(z.bbox.w,z.bbox.h));if(ss(v.bbox,z.bbox)>.28||q<O*.18){V=!0;break}}}if(!V){const U=P.length;if(P.push(v),A(v,U),P.length>=i)break}}return s==null||s("Detecting candidates: deduplicating...",1),P}function Hl(n,t,e){const i=new Uint8Array(n.length);for(let r=0;r<e;r++)for(let s=0;s<t;s++){let o=0,a=0;for(let l=-1;l<=1;l++){const u=r+l;if(!(u<0||u>=e))for(let h=-1;h<=1;h++){const c=s+h;c<0||c>=t||(o+=n[u*t+c],a++)}}i[r*t+s]=Math.round(o/Math.max(1,a))}return i}function ql(n,t,e,i){const r=Math.max(t,e);if(r<=i)return{gray:n,width:t,height:e,scale:1};const s=i/r,o=Math.max(1,Math.round(t*s)),a=Math.max(1,Math.round(e*s)),l=new Uint8Array(o*a);for(let u=0;u<a;u++){const h=Math.min(e-1,Math.floor(u/s));for(let c=0;c<o;c++){const y=Math.min(t-1,Math.floor(c/s));l[u*o+c]=n[h*t+y]}}return{gray:l,width:o,height:a,scale:s}}function ss(n,t){const e=Math.max(n.x,t.x),i=Math.max(n.y,t.y),r=Math.min(n.x+n.w,t.x+t.w),s=Math.min(n.y+n.h,t.y+t.h),o=Math.max(0,r-e),a=Math.max(0,s-i),l=o*a;if(l<=0)return 0;const u=n.w*n.h+t.w*t.h-l;return u>0?l/u:0}function as(n){const t=n.length;if(t===0)return{count:0,mean:0,std:1/0};let e=0;for(const s of n)e+=s;const i=e/t;let r=0;for(const s of n){const o=s-i;r+=o*o}return r/=t,{count:t,mean:i,std:Math.sqrt(Math.max(0,r))}}function Ql(n,t,e,i){return{p1:{x:n.x-t*i,y:n.y-e*i},p2:{x:n.x+t*i,y:n.y+e*i}}}function os(n,t,e,i,r){return[{x:n.p1.x+t*i,y:n.p1.y+e*i},{x:n.p2.x+t*i,y:n.p2.y+e*i},{x:n.p2.x+t*r,y:n.p2.y+e*r},{x:n.p1.x+t*r,y:n.p1.y+e*r}]}function Kl(n,t,e){let i=0;for(let r=0;r<4;r++){const s=e[r],o=e[(r+1)%4],a=(o.x-s.x)*(t-s.y)-(o.y-s.y)*(n-s.x);if(Math.abs(a)<=1e-6)continue;const l=a>0?1:-1;if(i===0)i=l;else if(i!==l)return!1}return!0}function ls(n,t,e,i){const r=xt(Nt(i,1),t,e);if(!r)return[];const s=[];for(let o=r.y;o<r.y+r.h;o++)for(let a=r.x;a<r.x+r.w;a++)Kl(a,o,i)&&s.push(n[o*t+a]);return s}function $l(n,t,e,i,r,s,o,a,l,u,h,c,y,f){const p=h*2,d=c*2,g=Math.hypot(i[1].x-i[0].x,i[1].y-i[0].y),x=Math.hypot(i[2].x-i[1].x,i[2].y-i[1].y),m=Math.hypot(i[2].x-i[3].x,i[2].y-i[3].y),b=Math.hypot(i[3].x-i[0].x,i[3].y-i[0].y),P=Math.max(...[g,x,m,b]),M=Math.max(2,Math.min(p,d)),w=Math.max(4,P*.25),C=Math.max(2,Math.min(12,M*.22)),F=Math.max(1,Math.min(C,Math.max(1,M*.5-1))),A=1,k=Math.max(8,Math.round(Math.min(w,C*3))),v=[[i[0],i[1],i[1],i[0]],[i[1],i[2],i[2],i[1]],[i[2],i[3],i[3],i[2]],[i[3],i[0],i[0],i[3]]],S=[[i[0],i[1],i[1],i[0]],[i[1],i[2],i[2],i[1]],[i[2],i[3],i[3],i[2]],[i[3],i[0],i[0],i[3]]],T=[],L=[],E=[{corners:[i[0],i[1]],seedLine:r,sideLength:g},{corners:[i[1],i[2]],seedLine:s,sideLength:x},{corners:[i[2],i[3]],seedLine:o,sideLength:m},{corners:[i[3],i[0]],seedLine:a,sideLength:b}];for(let G=0;G<E.length;G++){const X=E[G],st=Math.max(1,X.sideLength*.5-1),$=Math.max(1,Math.min(st,w*.5)),Z={x:(X.corners[0].x+X.corners[1].x)*.5,y:(X.corners[0].y+X.corners[1].y)*.5},lt=Ql(Z,X.seedLine.dirX,X.seedLine.dirY,$),ut=oe(n,t,e,lt.p1,lt.p2,$,Math.max(4,A+Math.max(C,F)+2)),nt=(ut==null?void 0:ut.line)||lt,Zt=nt.p2.x-nt.p1.x,Yt=nt.p2.y-nt.p1.y,Y=Math.hypot(Zt,Yt);if(!Number.isFinite(Y)||Y<=1e-6)return{ok:!1,score:0,meanSpread:1/0,meanSpreadLimit:1/0,avgStd:1/0,avgStdLimit:1/0,contrast:0,outerMean:0,outerSideMeans:[0,0,0,0],outerSideStds:[1/0,1/0,1/0,1/0],outerSideStdLimit:1/0,outerSideQuads:v,innerSideOk:!1,innerSideStds:[1/0,1/0,1/0,1/0],innerSideStdLimit:1/0,innerSideQuads:S};let pt=-Yt/Y,mt=Zt/Y;const gt={x:(nt.p1.x+nt.p2.x)*.5,y:(nt.p1.y+nt.p2.y)*.5};(gt.x-l)*pt+(gt.y-u)*mt<0&&(pt=-pt,mt=-mt);const ht=os(nt,pt,mt,A,A+C),fe=os(nt,pt,mt,-A,-(A+F));v[G]=ht,S[G]=fe,T.push(ls(n,t,e,ht)),L.push(ls(n,t,e,fe))}const V=T.map(as);if(V.some(G=>G.count<k||!Number.isFinite(G.std)))return{ok:!1,score:0,meanSpread:1/0,meanSpreadLimit:1/0,avgStd:1/0,avgStdLimit:1/0,contrast:0,outerMean:0,outerSideMeans:[0,0,0,0],outerSideStds:[1/0,1/0,1/0,1/0],outerSideStdLimit:1/0,outerSideQuads:v,innerSideOk:!1,innerSideStds:[1/0,1/0,1/0,1/0],innerSideStdLimit:1/0,innerSideQuads:S};const U=L.map(as);if(U.some(G=>G.count<k||!Number.isFinite(G.std)||!Number.isFinite(G.mean)))return{ok:!1,score:0,meanSpread:1/0,meanSpreadLimit:1/0,avgStd:1/0,avgStdLimit:1/0,contrast:0,outerMean:0,outerSideMeans:[0,0,0,0],outerSideStds:[1/0,1/0,1/0,1/0],outerSideStdLimit:1/0,outerSideQuads:v,innerSideOk:!1,innerSideStds:[1/0,1/0,1/0,1/0],innerSideStdLimit:1/0,innerSideQuads:S};const N=V.map(G=>G.mean),D=N.reduce((G,X)=>G+X,0)/N.length,B=U.reduce((G,X)=>G+X.mean,0)/U.length,z=Math.abs(B-D),q=Math.max(...N)-Math.min(...N),O=V.reduce((G,X)=>G+X.std,0)/V.length,W=Math.max(0,f),j=Math.max(6,Math.min(20,z*.45)),K=N,et=V.map(G=>G.std),tt=Math.max(j,Math.min(30,j*y)),ot=U.map(G=>G.std),rt=ot.every(G=>G<=tt),R=q<=W&&O<=j,Q=1/(1+q/Math.max(1,W)+O/Math.max(1,j));return{ok:R,score:Q,meanSpread:q,meanSpreadLimit:W,avgStd:O,avgStdLimit:j,contrast:z,outerMean:D,outerSideMeans:K,outerSideStds:et,outerSideStdLimit:j,outerSideQuads:v,innerSideOk:rt,innerSideStds:ot,innerSideStdLimit:tt,innerSideQuads:S}}function oe(n,t,e,i,r,s,o){const a=r.x-i.x,l=r.y-i.y,u=Math.hypot(a,l);if(!Number.isFinite(u)||u<=1e-6)return null;const h=a/u,c=l/u,y=-c,f=h,p=(i.x+r.x)*.5,d=(i.y+r.y)*.5,g=he({p1:i,p2:r},o+2);if(!xt(Nt(g||[i,r],2),t,e))return null;const m=Math.max(8,Math.round(s*2)+1),b=Math.max(8,Math.round(o*2)+1),_=m>1?s*2/(m-1):0,P=b>1?o*2/(b-1):0,M=Array.from({length:m},()=>new Array(b).fill(0));for(let v=0;v<m;v++){const S=-s+_*v;for(let T=0;T<b;T++){const L=-o+P*T,E=p+S*h+L*y,V=d+S*c+L*f;M[v][T]=Jl(n,t,e,E,V)}}const w=Gs(M,-o,-s,P,_,!0);if(w.length<8)return null;const C=w.map(v=>{const S=v.x,T=v.y;return{x:p+T*h+S*y,y:d+T*c+S*f,weight:v.weight}}),F=be(C);if(!F)return null;let A=F.dirX,k=F.dirY;return A*h+k*c<0&&(A=-A,k=-k),{line:{p1:{x:F.pointX-A*s,y:F.pointY-k*s},p2:{x:F.pointX+A*s,y:F.pointY+k*s}},fitPoints:C.map(v=>({x:v.x,y:v.y}))}}function Jl(n,t,e,i,r){if(t<=0||e<=0||n.length!==t*e)return 0;const s=Math.max(0,Math.min(t-1,i)),o=Math.max(0,Math.min(e-1,r)),a=Math.floor(s),l=Math.floor(o),u=Math.min(t-1,a+1),h=Math.min(e-1,l+1),c=s-a,y=o-l,f=n[l*t+a],p=n[l*t+u],d=n[h*t+a],g=n[h*t+u],x=f*(1-c)+p*c,m=d*(1-c)+g*c;return x*(1-y)+m*y}function Zl(n,t,e,i,r){if(i<=0||r<=0||i>=t-1||r>=e-1)return{gx:0,gy:0};const s=r*t+i;return{gx:(n[s+1]-n[s-1])*.5,gy:(n[s+t]-n[s-t])*.5}}function tc(n){if(n.length<20)return null;const t=n.map(w=>Math.max(0,w.weight));let e=0;for(const w of t)e=Math.max(e,w);if(!(e>0))return null;for(let w=0;w<t.length;w++)t[w]/=e;const i=w=>{let C=0,F=0,A=0;for(let j=0;j<n.length;j++){const K=w[j];K>0&&(C+=K,F+=n[j].x*K,A+=n[j].y*K)}if(!(C>0))return null;F/=C,A/=C;let k=0,v=0,S=0;for(let j=0;j<n.length;j++){const K=w[j];if(!(K>0))continue;const et=n[j].x-F,tt=n[j].y-A;k+=K*et*et,v+=K*et*tt,S+=K*tt*tt}k/=C,v/=C,S/=C;const T=k+S,L=k*S-v*v,E=-T,V=L,U=Math.max(0,E*E-4*V),N=-.5*(E+(E>=0?1:-1)*Math.sqrt(U)),D=Math.abs(N)>1e-12?N:0,B=Math.abs(N)>1e-12?V/N:T,z=Math.max(D,B);let q=0,O=1;Math.abs(v)>1e-10?(q=z-S,O=v):k>S&&(q=1,O=0);const W=Math.atan2(-q,O);return{centroid:{x:F,y:A},angle:W,totalWeight:C}},r=i(t);if(!r)return null;const s=Math.cos(r.angle),o=Math.sin(r.angle),a=new Array(2*16*8).fill(0);for(let w=0;w<n.length;w++){const C=n[w].x-r.centroid.x,F=n[w].y-r.centroid.y,A=C*s+F*o,k=Math.round(A*8+16*8);if(k>=3&&k<a.length-3)for(let v=-3;v<=3;v++)a[k+v]+=t[w]}let l=16*8;for(let w=-5*8+16*8;w<=5*8+16*8;w++)a[w]>a[l]&&(l=w);let u=l-1;for(;u>1&&a[u]>.05*a[l];)u--;let h=l+1;for(;h<a.length-1&&a[h]>.05*a[l];)h++;let c=Math.max(1,u-8);for(;c>1&&a[c]<=a[u];)c--;let y=Math.min(a.length-1,h+8);for(;y<a.length-1&&a[y]<=a[h];)y++;const f=a.slice();for(let w=1;w<f.length;w++)f[w]+=f[w-1];const p=f[f.length-1];if(!(p>0))return null;let d=0;for(let w=1;w<f.length;w++)Math.abs(f[w]-.1*p)<Math.abs(f[d]-.1*p)&&(d=w);let g=f.length-1;for(let w=f.length-2;w>0;w--)Math.abs(f[w]-.9*p)<Math.abs(f[g]-.9*p)&&(g=w);let x=d/8-16,m=g/8-16;const b=m-x;x-=b*.7,m+=b*.7,x=Math.max((c+u)/16-16,x),m=Math.min((y+h)/16-16,m);const _=t.slice();for(let w=0;w<n.length;w++){const C=n[w].x-r.centroid.x,F=n[w].y-r.centroid.y,A=C*s+F*o;_[w]=A>=x&&A<=m?t[w]**4*(1/(10+Math.abs(A))):0}const P=i(_);if(!P)return null;const M=[];for(let w=0;w<n.length;w++)_[w]>0&&M.push({x:n[w].x,y:n[w].y,weight:_[w]});return M.length<8?null:{centroid:P.centroid,angle:P.angle,keptSamples:M}}function ec(n,t,e,i,r,s=Ft){var C;const o=r.x-i.x,a=r.y-i.y,l=Math.hypot(o,a);if(!Number.isFinite(l)||l<=12)return null;const u=o/l,h=a/l,c=5,y=4*s+.5,f=(F,A,k,v,S)=>{const T={x:F.x-A*l*.5,y:F.y-k*l*.5},L={p1:T,p2:{x:T.x+A*l,y:T.y+k*l}},E=he(L,y+2),V=xt(Nt(E??[L.p1,L.p2],3),t,e),U=[],N=new Map;if(!V)return{reduced:null,scanlines:N};for(let D=V.y;D<V.y+V.h;D++)for(let B=V.x;B<V.x+V.w;B++){const z=B,q=D,O=z-T.x,W=q-T.y,j=O*A+W*k;if(!(j>c&&j<l-c))continue;const K=z-F.x,et=q-F.y,tt=K*v+et*S;if(Math.abs(tt)<12){const{gx:ot,gy:rt}=Zl(n,t,e,B,D),R=ot*ot+rt*rt;R>0&&U.push({x:z,y:q,weight:R})}if(Math.abs(tt)<y){const ot=N.get(D);ot?(B<ot.start&&(ot.start=B),B>ot.end&&(ot.end=B)):N.set(D,{start:B,end:B})}}return{reduced:tc(U),scanlines:N}};let p={x:(i.x+r.x)*.5,y:(i.y+r.y)*.5},d=u,g=h,x=-g,m=d,b=f(p,d,g,x,m);if(!b.reduced)return null;p=b.reduced.centroid,x=Math.cos(b.reduced.angle),m=Math.sin(b.reduced.angle),d=-m,g=x,d*u+g*h<0&&(d=-d,g=-g,x=-x,m=-m);let _=f(p,d,g,x,m);if(!_.reduced)return null;const P=Math.hypot(_.reduced.centroid.x-p.x,_.reduced.centroid.y-p.y);p=_.reduced.centroid,x=Math.cos(_.reduced.angle),m=Math.sin(_.reduced.angle),d=-m,g=x,d*u+g*h<0&&(d=-d,g=-g,x=-x,m=-m);let M=_;if(P>1){const F=f(p,d,g,x,m);F.reduced&&(M=F,p=F.reduced.centroid,x=Math.cos(F.reduced.angle),m=Math.sin(F.reduced.angle),d=-m,g=x,d*u+g*h<0&&(d=-d,g=-g))}const w=(((C=M.reduced)==null?void 0:C.keptSamples)??[]).map(F=>({x:F.x,y:F.y}));return w.length<8?null:{line:{p1:{x:p.x-d*l*.5,y:p.y-g*l*.5},p2:{x:p.x+d*l*.5,y:p.y+g*l*.5}},fitPoints:w,correctedScanlines:M.scanlines}}function nc(n,t,e){var b;const i=n.length,r=((b=n[0])==null?void 0:b.length)??0;if(r===0||i===0)return 0;const s=Math.max(0,Math.min(r-1,t)),o=Math.max(0,Math.min(i-1,e)),a=Math.floor(s),l=Math.floor(o),u=Math.min(r-1,a+1),h=Math.min(i-1,l+1),c=s-a,y=o-l,f=n[l][a],p=n[l][u],d=n[h][a],g=n[h][u],x=f*(1-c)+p*c,m=d*(1-c)+g*c;return x*(1-y)+m*y}function ic(n,t,e,i,r,s,o){var U;const a=i.p2.x-i.p1.x,l=i.p2.y-i.p1.y,u=Math.hypot(a,l);if(!Number.isFinite(u)||u<=1e-6)return null;const h=a/u,c=l/u,y=-c,f=h,p=(i.p1.x+i.p2.x)*.5,d=(i.p1.y+i.p2.y)*.5,g=he(i,s+2),x=xt(Nt(g||[i.p1,i.p2],3),t,e);if(!x)return null;const m=ll(n,t,e,x,0,0,cn(o==null?void 0:o.bayerPattern,"RAW edge refinement"),o==null?void 0:o.greenPhase,o==null?void 0:o.blackLevel),b=m.length;if((((U=m[0])==null?void 0:U.length)??0)<6||b<6)return null;const P=Math.max(8,Math.round(r*2)+1),M=Math.max(8,Math.round(s*2)+1),w=P>1?r*2/(P-1):0,C=M>1?s*2/(M-1):0,F=Array.from({length:P},()=>new Array(M).fill(0));for(let N=0;N<P;N++){const D=-r+w*N;for(let B=0;B<M;B++){const z=-s+C*B,q=p+D*h+z*y,O=d+D*c+z*f;F[N][B]=nc(m,q-x.x,O-x.y)}}const{gx:A,gy:k}=cl(F),v=A>=k,S=Gs(F,-s,-r,C,w,v);if(S.length<8)return null;const T=S.map(N=>{const D=N.x,B=N.y;return{x:p+B*h+D*y,y:d+B*c+D*f,weight:N.weight}}),L=be(T);if(!L)return null;let E=L.dirX,V=L.dirY;return E*h+V*c<0&&(E=-E,V=-V),{line:{p1:{x:L.pointX-E*r,y:L.pointY-V*r},p2:{x:L.pointX+E*r,y:L.pointY+V*r}},fitPoints:T.map(N=>({x:N.x,y:N.y}))}}function _e(n){const t=Math.max(0,Math.min(1,n));return t<=.04045?t/12.92:Math.pow((t+.055)/1.055,2.4)}function mr(n,t,e){if(t<=0||e<=0||n.length!==t*e)return new Uint8Array(Math.max(0,t*e));let i=1/0,r=-1/0;for(let p=0;p<n.length;p++){const d=n[p];Number.isFinite(d)&&(d<i&&(i=d),d>r&&(r=d))}if(!Number.isFinite(i)||!Number.isFinite(r)||r<=i+1e-9)return new Uint8Array(t*e);const s=1024,o=new Uint32Array(s),a=r-i;for(let p=0;p<n.length;p++){const d=Math.max(0,Math.min(1,(n[p]-i)/a)),g=Math.min(s-1,Math.max(0,Math.floor(d*(s-1))));o[g]++}const l=n.length,u=p=>{const d=l*p;let g=0;for(let x=0;x<s;x++)if(g+=o[x],g>=d)return i+x/Math.max(1,s-1)*a;return r},h=u(.01),c=u(.99),y=Math.max(1e-9,c-h),f=new Uint8Array(t*e);for(let p=0;p<n.length;p++){const d=Math.max(0,Math.min(1,(n[p]-h)/y));f[p]=Math.round(d*255)}return f}function Jn(n,t,e=0){const i=new Float32Array(n.width*n.height),r=n.data;for(let s=0,o=0;s<r.length;s+=4,o++)i[o]=Hs(r,s,t,e);return mr(i,n.width,n.height)}function rc(n){return Number.isFinite(n)?Math.max(0,Math.min(65535,Number(n))):0}function Hs(n,t,e,i=0){let r=n[t]/255,s=n[t+1]/255,o=n[t+2]/255;e&&(r=_e(r),s=_e(s),o=_e(o));const a=.2126*r+.7152*s+.0722*o;return Math.max(0,a-rc(i)/65535)}function sc(n,t){const e=n.width,i=n.height,r=n.data;if(r.length<e*i*3)return new Uint8Array(e*i);const s=new Float32Array(e*i);for(let o=0;o<e*i;o++){const a=o*3;t!==void 0?s[o]=r[a+t]:s[o]=.2126*r[a]+.7152*r[a+1]+.0722*r[a+2]}return mr(s,e,i)}function ac(n){const t=new Float32Array(n.width*n.height);for(let e=0;e<n.data.length;e++)t[e]=n.data[e];return mr(t,n.width,n.height)}function cs(n,t,e){const i=xt(t,n.width,n.height);if(!i)return null;const r=new Uint16Array(i.w*i.h*3),s=n.data;let o=0;for(let a=i.y;a<i.y+i.h;a++)for(let l=i.x;l<i.x+i.w;l++){const u=(a*n.width+l)*4;let h=s[u]/255,c=s[u+1]/255,y=s[u+2]/255;e&&(h=_e(h),c=_e(c),y=_e(y)),r[o++]=Math.max(0,Math.min(65535,Math.round(h*65535))),r[o++]=Math.max(0,Math.min(65535,Math.round(c*65535))),r[o++]=Math.max(0,Math.min(65535,Math.round(y*65535)))}return{data:r,width:i.w,height:i.h}}function us(n,t,e,i=0){const r=xt(t,n.width,n.height);if(!r)return null;const s=new Uint16Array(r.w*r.h*3),o=n.data;let a=0;for(let l=r.y;l<r.y+r.h;l++)for(let u=r.x;u<r.x+r.w;u++){const h=(l*n.width+u)*4,c=Math.max(0,Math.min(65535,Math.round(Hs(o,h,e,i)*65535)));s[a++]=c,s[a++]=c,s[a++]=c}return{data:s,width:r.w,height:r.h}}function oc(n,t){const e=xt(t,n.width,n.height);if(!e)return null;const i=new Uint16Array(e.w*e.h),r=n.data;let s=0;for(let o=e.y;o<e.y+e.h;o++)for(let a=e.x;a<e.x+e.w;a++){const l=(o*n.width+a)*4;i[s++]=Math.max(0,Math.min(65535,Math.round((.2126*r[l]+.7152*r[l+1]+.0722*r[l+2])*257)))}return{data:i,width:e.w,height:e.h}}function lc(n,t){const e=xt(t,n.width,n.height);if(!e)return null;const i=new Uint16Array(e.w*e.h);let r=0;for(let s=e.y;s<e.y+e.h;s++){const o=s*n.width;for(let a=e.x;a<e.x+e.w;a++)i[r++]=n.data[o+a]}return{data:i,width:e.w,height:e.h}}function Gt(n,t,e){return{x:n.x*t,y:n.y*e}}function cc(n,t,e){return{p1:Gt(n.p1,t,e),p2:Gt(n.p2,t,e)}}function Ki(n,t){const e=t(n);return{x:Number.isFinite(e.x)?e.x:n.x,y:Number.isFinite(e.y)?e.y:n.y}}function uc(n,t){return n.map(e=>Ki(e,t))}function ae(n,t,e,i=0,r=0){if(!n||n.length<8)return;const s=n.map(o=>({x:o.x*t-i,y:o.y*e-r})).filter(o=>Number.isFinite(o.x)&&Number.isFinite(o.y));return s.length>=8?s:void 0}function hc(n,t){return{p1:Ki(n.p1,t),p2:Ki(n.p2,t)}}function hn(n,t,e){return{p1:{x:n.p1.x-t,y:n.p1.y-e},p2:{x:n.p2.x-t,y:n.p2.y-e}}}function dc(n,t,e,i){const r=Math.max(0,Math.min(n.width-1,t)),o=(Math.max(0,Math.min(n.height-1,e))*n.width+r)*4;let a=n.data[o]/255,l=n.data[o+1]/255,u=n.data[o+2]/255;return i&&(a=_e(a),l=_e(l),u=_e(u)),(.2126*a+.7152*l+.0722*u)*65535}function qs(n){return n.kind==="u16-mono"}function en(n){return n.width}function nn(n){return n.height}function ci(n,t,e,i){if(qs(n)){const r=Math.max(0,Math.min(n.width-1,t)),s=Math.max(0,Math.min(n.height-1,e));return n.data[s*n.width+r]}return dc(n,t,e,i)}function fc(n,t,e,i){if(qs(n)&&n.coordinateSpace==="distorted-padded"){const r=Math.round(n.paddingOffsetX??0),s=Math.round(n.paddingOffsetY??0);return ci(n,t+r,e+s,i)}return ci(n,t,e,i)}function pc(n,t,e,i,r=3){const o=[...n,{x:(n[0].x+n[1].x+n[2].x+n[3].x)*.25,y:(n[0].y+n[1].y+n[2].y+n[3].y)*.25},{x:(n[0].x+n[1].x)*.5,y:(n[0].y+n[1].y)*.5},{x:(n[1].x+n[2].x)*.5,y:(n[1].y+n[2].y)*.5},{x:(n[2].x+n[3].x)*.5,y:(n[2].y+n[3].y)*.5},{x:(n[3].x+n[0].x)*.5,y:(n[3].y+n[0].y)*.5}].map(a=>we(a,t)).filter(a=>Number.isFinite(a.x)&&Number.isFinite(a.y));return o.length===0?null:xt(Nt(o,r),e,i)}function gr(n,t,e,i){const r=new Map;for(let s=n.y;s<n.y+n.h;s++)for(let o=n.x;o<n.x+n.w;o++){const a=Tt({x:o,y:s},t);if(!Number.isFinite(a.x)||!Number.isFinite(a.y))continue;const l=Math.round(a.x),u=Math.round(a.y);if(l<0||u<0||l>=e||u>=i)continue;const h=r.get(u);h?(l<h.start&&(h.start=l),l>h.end&&(h.end=l)):r.set(u,{start:l,end:l})}return r}function mc(n,t,e,i,r,s,o){const a=new Map,l=t.p2.x-t.p1.x,u=t.p2.y-t.p1.y,h=Math.hypot(l,u);if(!Number.isFinite(h)||h<=1e-6)return a;const c=l/h,y=u/h,f=-y,p=c,d={x:(t.p1.x+t.p2.x)*.5,y:(t.p1.y+t.p2.y)*.5},g=Math.max(1,e+1),x=Math.max(1,i+1.5);for(let m=n.y;m<n.y+n.h;m++)for(let b=n.x;b<n.x+n.w;b++){const _=b+.5,P=m+.5,M=_-d.x,w=P-d.y,C=M*c+w*y;if(!Number.isFinite(C)||Math.abs(C)>g)continue;const F=M*f+w*p;if(!Number.isFinite(F)||Math.abs(F)>x)continue;const A=Tt({x:_,y:P},r);if(!Number.isFinite(A.x)||!Number.isFinite(A.y))continue;const k=Math.round(A.x),v=Math.round(A.y);if(k<0||v<0||k>=s||v>=o)continue;const S=a.get(v);S?(k<S.start&&(S.start=k),k>S.end&&(S.end=k)):a.set(v,{start:k,end:k})}return a}function Qs(n,t,e,i,r,s){const o=new Map,a=t.p2.x-t.p1.x,l=t.p2.y-t.p1.y,u=Math.hypot(a,l);if(!Number.isFinite(u)||u<=1e-6)return o;const h=a/u,c=l/u,y=-c,f=h,p={x:(t.p1.x+t.p2.x)*.5,y:(t.p1.y+t.p2.y)*.5},d=Math.max(1,e+1),g=Math.max(1,i+1.5),x=xt(n,r,s);if(!x)return o;for(let m=x.y;m<x.y+x.h;m++)for(let b=x.x;b<x.x+x.w;b++){const _=b-p.x,P=m-p.y,M=_*h+P*c;if(!Number.isFinite(M)||Math.abs(M)>d)continue;const w=_*y+P*f;if(!Number.isFinite(w)||Math.abs(w)>g)continue;const C=o.get(m);C?(b<C.start&&(C.start=b),b>C.end&&(C.end=b)):o.set(m,{start:b,end:b})}return o}function Ks(n,t,e,i){const r=new Map;for(const[s,o]of n)for(let a=o.start;a<=o.end;a++){const l=we({x:a,y:s},t);if(!Number.isFinite(l.x)||!Number.isFinite(l.y))continue;const u=Math.round(l.x),h=Math.round(l.y);if(u<0||h<0||u>=e||h>=i)continue;const c=r.get(h);c?(u<c.start&&(c.start=u),u>c.end&&(c.end=u)):r.set(h,{start:u,end:u})}return r}function yr(n){return Math.abs(n.k1)<1e-4&&Math.abs(n.k2)<1e-4}function gc(n){return[{x:n.x,y:n.y},{x:n.x+n.w,y:n.y},{x:n.x+n.w,y:n.y+n.h},{x:n.x,y:n.y+n.h}]}function Xt(n,t,e,i,r){return we({x:i.x+n*t,y:i.y+n*e},r)}function xr(n,t,e,i,r){const o=Xt(n,t,e,i,r),a=Xt(n+1e-4,t,e,i,r);return{x:(a.x-o.x)/1e-4,y:(a.y-o.y)/1e-4}}function $s(n,t,e,i,r,s){let o=.01;const a=h=>{const c=Xt(h,t,e,i,s);return Math.hypot(c.x-r.x,c.y-r.y)},l=a(n),u=a(n+o);if(!Number.isFinite(l)||!Number.isFinite(u))return null;if(l>u){let h=n,c=n+o;for(let y=0;y<24;y++){o*=2;const f=h+o,p=a(f),d=a(c);if(!Number.isFinite(p)||!Number.isFinite(d))break;if(p>=d)return{a:h,b:f};h=c,c=f}}else{let h=n,c=n+o;for(let y=0;y<24;y++){o*=2;const f=c-o,p=a(f),d=a(h);if(!Number.isFinite(p)||!Number.isFinite(d))break;if(p>=d)return{a:f,b:c};c=h,h=f}}return{a:n-Math.max(.5,o),b:n+Math.max(.5,o)}}function yc(n,t,e=33){const i=n.p2.x-n.p1.x,r=n.p2.y-n.p1.y,s=Math.hypot(i,r);if(!Number.isFinite(s)||s<=1e-6)return[we(n.p1,t),we(n.p2,t)];const o=i/s,a=r/s,l={x:(n.p1.x+n.p2.x)*.5,y:(n.p1.y+n.p2.y)*.5},u=s*.5,h=Math.max(9,e),c=[];for(let y=0;y<h;y++){const f=h===1?.5:y/(h-1),p=-u+f*(u*2);c.push(Xt(p,o,a,l,t))}return c}function xc(n,t,e,i,r,s,o=1){const a=n.p2.x-n.p1.x,l=n.p2.y-n.p1.y,u=Math.hypot(a,l);if(!Number.isFinite(u)||u<=1e-6)return null;const h=a/u,c=l/u,y={x:(n.p1.x+n.p2.x)*.5,y:(n.p1.y+n.p2.y)*.5},f=Math.max(24,Math.round(e*2)+1),p=[];for(let d=0;d<f;d++){const g=f===1?.5:d/(f-1),x=-e+g*(e*2),m=Xt(x,h,c,y,t),b=xr(x,h,c,y,t),_=Math.hypot(b.x,b.y);if(!Number.isFinite(_)||_<=1e-9)continue;const P=-b.y/_,M=b.x/_;p.push({x:m.x+P*(i+o),y:m.y+M*(i+o)},{x:m.x-P*(i+o),y:m.y-M*(i+o)})}if(p.length<2){const d={p1:we(n.p1,t),p2:we(n.p2,t)},g=he(d,i+o);return g?xt(Nt(g,o),r,s):null}return xt(Nt(p,o),r,s)}function tn(n,t,e,i,r){return Tt({x:i.x+n*t,y:i.y+n*e},r)}function bc(n,t,e,i,r,s){let o=.01;const a=h=>{const c=tn(h,t,e,i,s);return Math.hypot(c.x-r.x,c.y-r.y)},l=a(n),u=a(n+o);if(!Number.isFinite(l)||!Number.isFinite(u))return null;if(l>u){let h=n,c=n+o;for(let y=0;y<24;y++){o*=2;const f=h+o,p=a(f),d=a(c);if(!Number.isFinite(p)||!Number.isFinite(d))break;if(p>=d)return{a:h,b:f};h=c,c=f}}else{let h=n,c=n+o;for(let y=0;y<24;y++){o*=2;const f=c-o,p=a(f),d=a(h);if(!Number.isFinite(p)||!Number.isFinite(d))break;if(p>=d)return{a:f,b:c};c=h,h=f}}return{a:n-Math.max(.5,o),b:n+Math.max(.5,o)}}function br(n,t,e){const i=(t.x-n.x)*(t.y-e.y)-(t.x-e.x)*(t.y-n.y);if(Math.abs(i)<=1e-12)return .5*(n.x+e.x);const r=(t.x-n.x)*(t.x-n.x)*(t.y-e.y)-(t.x-e.x)*(t.x-e.x)*(t.y-n.y),s=t.x-.5*r/i;return Number.isFinite(s)?s:.5*(n.x+e.x)}function _c(n,t,e,i,r){const o=tn(n,t,e,i,r),a=tn(n+1e-4,t,e,i,r);return{x:(a.x-o.x)/1e-4,y:(a.y-o.y)/1e-4}}function wc(n,t,e,i,r,s,o=!1,a){const l=[],u=[],h=a?Ft*2:Ft,c=Math.max(1,Math.min(s,h)),y=i.p2.x-i.p1.x,f=i.p2.y-i.p1.y,p=Math.hypot(y,f);if(!Number.isFinite(p)||p<=1e-6)return null;const d=y/p,g=f/p,x={x:(i.p1.x+i.p2.x)*.5,y:(i.p1.y+i.p2.y)*.5},m={p1:Tt(i.p1,e),p2:Tt(i.p2,e)},b=m.p2.x-m.p1.x,_=m.p2.y-m.p1.y,P=Math.hypot(b,_);if(!Number.isFinite(P)||P<=1e-6)return null;const M=b/P,w=_/P,C=-w,F=M,A={x:(m.p1.x+m.p2.x)*.5,y:(m.p1.y+m.p2.y)*.5},k=xt(a||Nt(he(i,s+2)??[i.p1,i.p2],2),n.width,n.height);if(!k)return null;const v=mc(k,i,r,c,e,en(t),nn(t));if(v.size===0)return null;const S=!yr(e);for(const[N,D]of v)if(!(N<0||N>=nn(t)))for(let B=D.start;B<=D.end;B++){if(B<0||B>=en(t))continue;const z={x:B,y:N};let q,O;if(S){const W=we(z,e);if(!Number.isFinite(W.x)||!Number.isFinite(W.y)||Math.round(W.x)<0||Math.round(W.x)>=n.width||Math.round(W.y)<0||Math.round(W.y)>=n.height)continue;const j=W.x-x.x,K=W.y-x.y,et=j*d+K*g;if(!Number.isFinite(et))continue;q=et,O=j*-g+K*d;const tt=bc(et,d,g,x,z,e);if(!tt)continue;const ot=.5*(tt.a+tt.b),rt=tn(tt.a,d,g,x,e),R=tn(ot,d,g,x,e),Q=tn(tt.b,d,g,x,e),G=br({x:tt.a,y:Math.hypot(rt.x-z.x,rt.y-z.y)},{x:ot,y:Math.hypot(R.x-z.x,R.y-z.y)},{x:tt.b,y:Math.hypot(Q.x-z.x,Q.y-z.y)});if(!Number.isFinite(G))continue;q=G;const X=_c(G,d,g,x,e),st=Math.hypot(X.x,X.y);if(!Number.isFinite(st)||st<=1e-9)continue;const $=X.x/st,lt=-(X.y/st),ut=$,nt=tn(G,d,g,x,e);O=(z.x-nt.x)*lt+(z.y-nt.y)*ut}else{const W=z.x-A.x,j=z.y-A.y;q=W*M+j*w,O=W*C+j*F}!Number.isFinite(q)||Math.abs(q)>r||!Number.isFinite(O)||Math.abs(O)>c||(l.push(O),u.push(ci(t,B,N,o)))}if(l.length<8)return null;const T=Tt(i.p1,e),L=Tt(i.p2,e),E=L.x-T.x,V=L.y-T.y,U=Math.abs(E)>=Math.abs(V)?1:2;return _n(l,u,U,h)}function Mc(n,t,e,i,r,s,o=!1,a,l,u,h=!1){const c=[],y=[],f=a?Ft*2:Ft,p=Math.max(1,Math.min(s,f)),d=i.p2.x-i.p1.x,g=i.p2.y-i.p1.y,x=Math.hypot(d,g);if(!Number.isFinite(x)||x<=1e-6)return null;const m=d/x,b=g/x,_=-b,P=m,M={x:(i.p1.x+i.p2.x)*.5,y:(i.p1.y+i.p2.y)*.5},w=xt(a||Nt(he(i,f*4+2)??[i.p1,i.p2],2),n.width,n.height);if(!w)return null;const C=l??(u?gr(xt(u,en(t),nn(t))??u,e,n.width,n.height):Qs(w,i,Math.max(1,r),p*4+.5,n.width,n.height));if(C.size===0)return null;const F=Ks(C,e,en(t),nn(t));if(F.size===0)return null;const A=!yr(e);for(const[v,S]of F)for(let T=S.start;T<=S.end;T++){const L={x:T,y:v},E=Tt(L,e);if(!Number.isFinite(E.x)||!Number.isFinite(E.y)||Math.round(E.x)<0||Math.round(E.x)>=n.width||Math.round(E.y)<0||Math.round(E.y)>=n.height)continue;const V=E.x-M.x,U=E.y-M.y,N=V*m+U*b;let D=V*_+U*P;if(A){const B=$s(N,m,b,M,L,e);if(!B)continue;const z=.5*(B.a+B.b),q=Xt(B.a,m,b,M,e),O=Xt(z,m,b,M,e),W=Xt(B.b,m,b,M,e),j=br({x:B.a,y:Math.hypot(q.x-L.x,q.y-L.y)},{x:z,y:Math.hypot(O.x-L.x,O.y-L.y)},{x:B.b,y:Math.hypot(W.x-L.x,W.y-L.y)});if(!Number.isFinite(j))continue;const K=xr(j,m,b,M,e),et=Math.hypot(K.x,K.y);if(!Number.isFinite(et)||et<=1e-9)continue;const tt=K.x/et,rt=-(K.y/et),R=tt,Q=Xt(j,m,b,M,e);D=(L.x-Q.x)*rt+(L.y-Q.y)*R}!Number.isFinite(N)||Math.abs(N)>Math.max(1,r)||!Number.isFinite(D)||Math.abs(D)>p||(c.push(D),y.push(fc(t,T,v,o)))}if(c.length<8)return null;const k=Math.abs(d)>=Math.abs(g)?1:2;return h?wn(c,y,k,f):_n(c,y,k,f)}function Sc(n,t,e,i,r,s){const o=n.width,a=n.height,l=cn(n.bayerPattern,"corrected RAW edge measurement"),u=s!=null&&s.correctedRect?Ft*2:Ft,h=Math.max(1,Math.min(r,u)),c=(s==null?void 0:s.restrictToStrip)??!0,y=e.p2.x-e.p1.x,f=e.p2.y-e.p1.y,p=Math.hypot(y,f);if(!Number.isFinite(p)||p<=1e-6)return null;const d=y/p,g=f/p,x=-g,m=d,b={x:(e.p1.x+e.p2.x)*.5,y:(e.p1.y+e.p2.y)*.5},_={p1:{x:b.x-d*Math.max(1,i),y:b.y-g*Math.max(1,i)},p2:{x:b.x+d*Math.max(1,i),y:b.y+g*Math.max(1,i)}},P=he(_,h+2),M=(s!=null&&s.fixedRawRect?xt(s.fixedRawRect,o,a):null)??(s!=null&&s.correctedRect?fr(s.correctedRect,t,o,a):null)??(P?pc(P,t,o,a,2):null);if(!M)return null;const w=[],C=[];for(let k=M.y;k<M.y+M.h;k++){const v=k*o;for(let S=M.x;S<M.x+M.w;S++){if(!_t(S,k,l,s==null?void 0:s.greenPhase))continue;const T=Tt({x:S,y:k},t);if(!Number.isFinite(T.x)||!Number.isFinite(T.y))continue;const L=T.x-b.x,E=T.y-b.y,V=L*d+E*g;if(!Number.isFinite(V)||c&&Math.abs(V)>Math.max(1,i))continue;const U=L*x+E*m;if(!Number.isFinite(U)||c&&Math.abs(U)>h)continue;w.push(U);let N;N=Math.max(0,n.data[v+S]-bn(s==null?void 0:s.blackLevel,S,k)),C.push(N)}}if(w.length<8)return null;const F=Math.abs(y)>=Math.abs(f)?1:2,A=Math.max(2,(s==null?void 0:s.shortSidePxOverride)??(c?h*2:Math.min(M.w,M.h)));return gi(w,C,A,s==null?void 0:s.manualBinSize,F,s==null?void 0:s.preferAutoPerEdgeBin,!1,!!(s!=null&&s.forceLegacyModel))}function Js(n,t,e,i,r,s=!1,o,a,l,u=!1){if(a&&l){const C=l.p2.x-l.p1.x,F=l.p2.y-l.p1.y,A=Math.hypot(C,F);if(Number.isFinite(A)&&A>1e-6)return wc(a,n,t,l,Math.max(1,A*.5),r,s,o)}const h=[],c=[],y=Ft,f=e.p2.x-e.p1.x,p=e.p2.y-e.p1.y,d=Math.hypot(f,p);if(!Number.isFinite(d)||d<=1e-6)return null;const g=f/d,x=p/d,m={x:(e.p1.x+e.p2.x)*.5,y:(e.p1.y+e.p2.y)*.5},b=(o?xt(o,en(n),nn(n)):null)??xc(e,t,i+1,r+1,en(n),nn(n),1);if(!b)return null;const _=gr(b,t,en(n),nn(n));if(_.size===0)return null;const P=-x,M=g;for(const[C,F]of _)for(let A=F.start;A<=F.end;A++){const k={x:A,y:C};let v=(k.x-m.x)*g+(k.y-m.y)*x,S=(k.x-m.x)*P+(k.y-m.y)*M;!Number.isFinite(v)||Math.abs(v)>i+1||!Number.isFinite(S)||Math.abs(S)>=y||(h.push(S),c.push(ci(n,A,C,s)))}if(h.length<8)return null;const w=Math.abs(f)>=Math.abs(p)?1:2;return u?wn(h,c,w,y):_n(h,c,w,y)}function vc(n,t,e){var y;const i=e.sourceMode??(t.isThreePlane?"three-plane":"rggb-raw"),r=e.useQuadraticProjection!==!1,s=!!e.forceRenderedMeasurement,o=n.width,a=n.height,l=e.threePlaneChannel,u=pr(e.detectionTuning),h=e.monochromeBlackLevel??0;if(i==="rggb-raw"&&!s){if(!t||t.isThreePlane)return null;const f=Nl(t,e.greenPhase),p=o/Math.max(1,t.width),d=a/Math.max(1,t.height);return{sourceMode:i,detectionGray:f,detectionWidth:t.width,detectionHeight:t.height,detectToMeasureX:1,detectToMeasureY:1,detectToDisplayX:p,detectToDisplayY:d,measureToDisplayX:p,measureToDisplayY:d,detectPointToDisplay:g=>Gt(g,p,d),measurePointToDisplay:g=>Gt(g,p,d),displayPointToDetect:g=>Gt(g,1/Math.max(1e-9,p),1/Math.max(1e-9,d)),measureUsesDisplayLine:!1,measureWidth:t.width,measureHeight:t.height,refineLine:(g,x,m)=>ic(t.data,t.width,t.height,{p1:g,p2:x},m*.5,Math.max(4,m*.2),{greenPhase:e.greenPhase,bayerPattern:t.bayerPattern})||oe(f,t.width,t.height,g,x,m*.5,Math.max(4,m*.2)),measureEdge:(g,x,m,b,_)=>xe(t.data,t.width,t.height,g,x,m,b,{greenOnly:!0,greenPhase:e.greenPhase,bayerPattern:t.bayerPattern,blackLevel:e.blackLevel??void 0,preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?ae(_==null?void 0:_.fitPoints,1,1):void 0})}}if(s){const f=!!e.distortionCurveApplied&&!!e.distortionModel,p=i==="rggb-raw"&&!!e.distortionCorrected&&!!e.distortionModel&&!t.isThreePlane,d=!!e.distortionCorrected&&!!e.distortionModel&&!!e.distortionOriginalSamplingPlane,g=!!e.distortionCorrected&&!!e.distortionSamplingPlane,x=n,m=Jn(x,!!e.sfrHasGamma,i==="unmix-bw"?h:0);return{sourceMode:i,detectionGray:m,detectionWidth:x.width,detectionHeight:x.height,detectToMeasureX:1,detectToMeasureY:1,detectToDisplayX:1,detectToDisplayY:1,measureToDisplayX:1,measureToDisplayY:1,detectPointToDisplay:b=>b,measurePointToDisplay:b=>b,displayPointToDetect:b=>b,measureUsesDisplayLine:!1,measureWidth:x.width,measureHeight:x.height,refineLine:(b,_,P)=>(d?ec(m,n.width,n.height,b,_,Ft):null)||oe(m,n.width,n.height,b,_,P*.5,Math.max(4,P*.2)),measureEdge:(b,_,P,M,w)=>{const C=e.distortionModel?fr(b,e.distortionModel,t.width,t.height):null;if(p){const v={p1:Tt(_.p1,e.distortionModel),p2:Tt(_.p2,e.distortionModel)},S=Math.hypot(v.p2.x-v.p1.x,v.p2.y-v.p1.y),T=Math.max(2,S*.5*u.sampleHalfWidthRatio);return Wo(t,e.distortionModel,v,Math.max(1,S*.5),T,{greenPhase:e.greenPhase,blackLevel:e.blackLevel??void 0,correctedRect:b})}if(p)return Sc(t,e.distortionModel,_,P,M,{greenPhase:e.greenPhase,blackLevel:e.blackLevel??void 0,correctedRect:b,fixedRawRect:C,preferAutoPerEdgeBin:!0});if(d)return Mc(n,e.distortionOriginalSamplingPlane,e.distortionModel,_,P,M,!!e.sfrHasGamma,b,(w==null?void 0:w.correctedScanlines)??null,C);if(f){const v={p1:Tt(_.p1,e.distortionModel),p2:Tt(_.p2,e.distortionModel)},S=Math.hypot(v.p2.x-v.p1.x,v.p2.y-v.p1.y);return Js(e.distortionOriginalSamplingPlane??e.distortionSamplingPlane??e.distortionSamplingImage??n,e.distortionModel,v,Math.max(1,S*.5),M,!!e.sfrHasGamma,b,e.distortionBaseImage??n,_)}if(g){const v=lc(e.distortionSamplingPlane,b);if(!v)return null;const S=hn(_,b.x,b.y);return xe(v.data,v.width,v.height,{x:0,y:0,w:v.width,h:v.height},S,P,M,{preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?ae(w==null?void 0:w.fitPoints,1,1,b.x,b.y):void 0})}const A=i==="unmix-bw"?us(n,b,!!e.sfrHasGamma,h):cs(n,b,!!e.sfrHasGamma);if(!A)return null;const k=hn(_,b.x,b.y);return xe(A.data,A.width,A.height,{x:0,y:0,w:A.width,h:A.height},k,P,M,{isThreePlane:!0,threePlaneChannel:void 0,preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?ae(w==null?void 0:w.fitPoints,1,1,b.x,b.y):void 0})}}}if(i==="three-plane"){if(t.isThreePlane&&!e.sfrHasGamma){const p=sc(t,l),d=o/Math.max(1,t.width),g=a/Math.max(1,t.height);return{sourceMode:i,detectionGray:p,detectionWidth:t.width,detectionHeight:t.height,detectToMeasureX:1,detectToMeasureY:1,detectToDisplayX:d,detectToDisplayY:g,measureToDisplayX:d,measureToDisplayY:g,detectPointToDisplay:x=>Gt(x,d,g),measurePointToDisplay:x=>Gt(x,d,g),displayPointToDetect:x=>Gt(x,1/Math.max(1e-9,d),1/Math.max(1e-9,g)),measureUsesDisplayLine:!1,measureWidth:t.width,measureHeight:t.height,refineLine:(x,m,b)=>oe(p,t.width,t.height,x,m,b*.5,Math.max(4,b*.2)),measureEdge:(x,m,b,_,P)=>xe(t.data,t.width,t.height,x,m,b,_,{isThreePlane:!0,threePlaneChannel:l,preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?ae(P==null?void 0:P.fitPoints,1,1):void 0})}}const f=Jn(n,!!e.sfrHasGamma);return{sourceMode:i,detectionGray:f,detectionWidth:n.width,detectionHeight:n.height,detectToMeasureX:1,detectToMeasureY:1,detectToDisplayX:1,detectToDisplayY:1,measureToDisplayX:1,measureToDisplayY:1,detectPointToDisplay:p=>p,measurePointToDisplay:p=>p,displayPointToDetect:p=>p,measureUsesDisplayLine:!1,measureWidth:n.width,measureHeight:n.height,refineLine:(p,d,g)=>oe(f,n.width,n.height,p,d,g*.5,Math.max(4,g*.2)),measureEdge:(p,d,g,x,m)=>{const b=cs(n,p,!!e.sfrHasGamma);if(!b)return null;const _=hn(d,p.x,p.y);return xe(b.data,b.width,b.height,{x:0,y:0,w:b.width,h:b.height},_,g,x,{isThreePlane:!0,threePlaneChannel:l,preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?ae(m==null?void 0:m.fitPoints,1,1,p.x,p.y):void 0})}}}if(i==="unmix-bw"){if(t&&!t.isThreePlane&&e.displaySettings){const p=so(t,e.displaySettings,e.blackLevel??e.monochromeBlackLevel??void 0);if(p){const d=ac(p),g=o/Math.max(1,t.width),x=a/Math.max(1,t.height);return{sourceMode:i,detectionGray:d,detectionWidth:t.width,detectionHeight:t.height,detectToMeasureX:1,detectToMeasureY:1,detectToDisplayX:g,detectToDisplayY:x,measureToDisplayX:g,measureToDisplayY:x,detectPointToDisplay:m=>Gt(m,g,x),measurePointToDisplay:m=>Gt(m,g,x),displayPointToDetect:m=>Gt(m,1/Math.max(1e-9,g),1/Math.max(1e-9,x)),measureUsesDisplayLine:!1,measureWidth:t.width,measureHeight:t.height,refineLine:(m,b,_)=>oe(d,t.width,t.height,m,b,_*.5,Math.max(4,_*.2)),measureEdge:(m,b,_,P,M)=>xe(p.data,p.width,p.height,m,b,_,P,{preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?ae(M==null?void 0:M.fitPoints,1,1):void 0})}}}const f=Jn(n,!!e.sfrHasGamma,h);return{sourceMode:i,detectionGray:f,detectionWidth:n.width,detectionHeight:n.height,detectToMeasureX:1,detectToMeasureY:1,detectToDisplayX:1,detectToDisplayY:1,measureToDisplayX:1,measureToDisplayY:1,detectPointToDisplay:p=>p,measurePointToDisplay:p=>p,displayPointToDetect:p=>p,measureUsesDisplayLine:!1,measureWidth:n.width,measureHeight:n.height,refineLine:(p,d,g)=>oe(f,n.width,n.height,p,d,g*.5,Math.max(4,g*.2)),measureEdge:(p,d,g,x,m)=>{const b=us(n,p,!!e.sfrHasGamma,h);if(!b)return null;const _=hn(d,p.x,p.y);return xe(b.data,b.width,b.height,{x:0,y:0,w:b.width,h:b.height},_,g,x,{isThreePlane:!0,threePlaneChannel:void 0,preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?ae(m==null?void 0:m.fitPoints,1,1,p.x,p.y):void 0})}}}const c=Jn(n,!1);if(t&&!t.isThreePlane&&((y=e.displaySettings)==null?void 0:y.renderMode)==="advanced-zero-dep"&&e.displaySettings.advancedZeroDep){const f=o/Math.max(1,t.width),p=a/Math.max(1,t.height);return{sourceMode:i,detectionGray:c,detectionWidth:n.width,detectionHeight:n.height,detectToMeasureX:t.width/Math.max(1,n.width),detectToMeasureY:t.height/Math.max(1,n.height),detectToDisplayX:1,detectToDisplayY:1,measureToDisplayX:f,measureToDisplayY:p,detectPointToDisplay:d=>d,measurePointToDisplay:d=>Gt(d,f,p),displayPointToDetect:d=>d,measureUsesDisplayLine:!1,measureWidth:t.width,measureHeight:t.height,refineLine:(d,g,x)=>oe(c,n.width,n.height,d,g,x*.5,Math.max(4,x*.2)),measureEdge:(d,g,x,m,b)=>{const _=ro(t,d,e.displaySettings);if(!_||_.width<8||_.height<8)return null;const P=hn(g,d.x,d.y);return xe(_.data,_.width,_.height,{x:0,y:0,w:_.width,h:_.height},P,x,m,{preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?ae(b==null?void 0:b.fitPoints,t.width/Math.max(1,n.width),t.height/Math.max(1,n.height),d.x,d.y):void 0})}}}if(t&&!t.isThreePlane){const f=o/Math.max(1,t.width),p=a/Math.max(1,t.height);return{sourceMode:i,detectionGray:c,detectionWidth:n.width,detectionHeight:n.height,detectToMeasureX:t.width/Math.max(1,n.width),detectToMeasureY:t.height/Math.max(1,n.height),detectToDisplayX:1,detectToDisplayY:1,measureToDisplayX:f,measureToDisplayY:p,detectPointToDisplay:d=>d,measurePointToDisplay:d=>Gt(d,f,p),displayPointToDetect:d=>d,measureUsesDisplayLine:!1,measureWidth:t.width,measureHeight:t.height,refineLine:(d,g,x)=>oe(c,n.width,n.height,d,g,x*.5,Math.max(4,x*.2)),measureEdge:(d,g,x,m,b)=>xe(t.data,t.width,t.height,d,g,x,m,{blackLevel:e.blackLevel??void 0,preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?ae(b==null?void 0:b.fitPoints,t.width/Math.max(1,n.width),t.height/Math.max(1,n.height)):void 0})}}return{sourceMode:i,detectionGray:c,detectionWidth:n.width,detectionHeight:n.height,detectToMeasureX:1,detectToMeasureY:1,detectToDisplayX:1,detectToDisplayY:1,measureToDisplayX:1,measureToDisplayY:1,detectPointToDisplay:f=>f,measurePointToDisplay:f=>f,displayPointToDetect:f=>f,measureUsesDisplayLine:!1,measureWidth:n.width,measureHeight:n.height,refineLine:(f,p,d)=>oe(c,n.width,n.height,f,p,d*.5,Math.max(4,d*.2)),measureEdge:(f,p,d,g,x)=>{const m=oc(n,f);if(!m)return null;const b=hn(p,f.x,f.y);return wl(m.data,m.width,m.height,{x:0,y:0,w:m.width,h:m.height},b,d,g,{blackLevel:e.monochromeBlackLevel??void 0,preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?ae(x==null?void 0:x.fitPoints,1,1,f.x,f.y):void 0})}}}function Pc(n,t,e){var u,h,c,y,f,p,d,g;if(!n||!t)return[];(u=e.onProgress)==null||u.call(e,"Preparing source context...",0);const i=vc(n,t,e);if(!i)return[];(h=e.onProgress)==null||h.call(e,"Preparing source context...",.08);const r=pr(e.detectionTuning),s=Math.min(1e3,Math.max(1,e.maxRegions??1e3)),o=Math.max(4,e.maxEdges??s*4);(c=e.onProgress)==null||c.call(e,"Detecting candidates...",.12);const a=jl(i.detectionGray,i.detectionWidth,i.detectionHeight,s,e.detectionTuning,(x,m)=>{var b;(b=e.onProgress)==null||b.call(e,x,.12+.08*Math.max(0,Math.min(1,m)))});(y=e.onProgress)==null||y.call(e,"Detecting candidates...",.2);const l=[];for(let x=0;x<a.length;x++){const m=a.length<=0?1:x/a.length;if((f=e.onProgress)==null||f.call(e,`Measuring edges: region ${x+1}/${a.length}`,.2+.72*Math.min(1,m)),l.length>=o)break;const b=a[x],_=b.corners,P=`auto-region-${x+1}`;for(let M=0;M<4&&((p=e.onProgress)==null||p.call(e,`Measuring edges: region ${x+1}/${a.length}, edge ${M+1}/4`,.2+.72*Math.min(1,(x+M/4)/Math.max(1,a.length))),!(l.length>=o));M=M+1){const w=_[M],C=_[(M+1)%4],F=C.x-w.x,A=C.y-w.y,k=Math.hypot(F,A);if(!Number.isFinite(k)||k<24)continue;const v=.125,S={x:w.x+F*v,y:w.y+A*v},T={x:C.x-F*v,y:C.y-A*v},L=Math.hypot(T.x-S.x,T.y-S.y);if(!Number.isFinite(L)||L<12)continue;const E=i.refineLine(S,T,L),V=(E!=null&&E.fitPoints?es(E.fitPoints):null)||(E==null?void 0:E.line)||{p1:S,p2:T},U=cc(V,i.detectToMeasureX,i.detectToMeasureY),N=uc((E==null?void 0:E.fitPoints)??[],i.detectPointToDisplay),D=(N.length>=2?es(N):null)||hc(U,i.measurePointToDisplay),B=i.measureUsesDisplayLine?D:U,z=B.p2.x-B.p1.x,q=B.p2.y-B.p1.y,O=Math.hypot(z,q);if(!Number.isFinite(O)||O<=1e-6)continue;const W=D.p2.x-D.p1.x,j=D.p2.y-D.p1.y,K=Math.hypot(W,j);if(!Number.isFinite(K)||K<=1e-6)continue;const et=!!e.distortionCurveApplied&&!!e.distortionModel,tt=W/K,ot=j/K;let rt=ot,R=-tt;const Q=(D.p1.x+D.p2.x)*.5,G=(D.p1.y+D.p2.y)*.5,X=i.detectPointToDisplay({x:b.centerX,y:b.centerY}),st=X.x,$=X.y;(Q-st)*rt+(G-$)*R<0&&(rt=-rt,R=-R);const Z=O*.5,lt=Math.max(2,O*r.sampleHalfWidthRatio),ut=Math.max(2,K*r.sampleHalfWidthRatio),nt=et?{p1:Tt(D.p1,e.distortionModel),p2:Tt(D.p2,e.distortionModel)}:void 0,Zt=nt?Math.max(1,Math.hypot(nt.p2.x-nt.p1.x,nt.p2.y-nt.p1.y)*.5):Z,Yt=et?D:U,Y=et?ut:lt,pt=he(Yt,Y);if(!pt)continue;const mt=xt(Nt(pt,2),et?n.width:i.measureWidth,et?n.height:i.measureHeight);if(!mt)continue;const gt=e.distortionCorrected&&e.distortionModel&&i.sourceMode==="rggb-raw"?fr(mt,e.distortionModel,t.width,t.height):null,ht=et?Js(e.distortionSamplingPlane??e.distortionSamplingImage??n,e.distortionModel,nt,Zt,Y,!!e.sfrHasGamma,mt,e.distortionBaseImage??n,Yt):i.measureEdge(mt,B,Z,Y,E);if(!ht||(ht.autoLikeUsed=!0,!bl(ht,e.useDeshading,0)))continue;const fe=e.useNR?-1:12,Wt=Al([ht],fe,null,e.useDeshading,0,!0);if(!Wt||Wt.mtf50===null||!_l(Wt.lsfCropped))continue;const Ot=nt?yc(nt,e.distortionModel,Math.max(21,Math.round(K*.5))):ht.quadraticProjectionUsed?gl(N,D,Math.max(21,Math.round(K*.5))):void 0,Re=nt&&Ot&&Ot.length>=2?Nt(Ot,ut+2):null,Le=Re?gc(Re):he(D,ut);if(!Le)continue;const Mn=e.distortionCorrected?mt:Re??Nt(Le,2);let un={x:Q+rt*(ut+12),y:G+R*(ut+12)},Rn=ns(tt,ot);if(Ot&&Ot.length>=3){const Sn=Math.floor(Ot.length/2),Ee=Ot[Math.max(0,Sn-1)],Ue=Ot[Math.min(Ot.length-1,Sn+1)],pe=Ot[Sn],Me=Ue.x-Ee.x,vn=Ue.y-Ee.y,De=Math.hypot(Me,vn);if(De>1e-6){const Be=vn/De,Oe=-Me/De;Rn=ns(Me/De,vn/De);const ze={x:pe.x-st,y:pe.y-$},Ve=ze.x*Be+ze.y*Oe>=0?1:-1;un={x:pe.x+Be*Ve*(ut+12),y:pe.y+Oe*Ve*(ut+12)}}}l.push({id:`${P}-edge-${M+1}`,regionId:P,sourceMode:i.sourceMode,edgeIndex:M,label:Wt.mtf50.toFixed(3),mtf50:Wt.mtf50,angle:Rn,orientation:ht.orientation,edgeData:ht,sourceRect:Mn,rawSourceRect:(i.sourceMode==="rggb-raw"?gt??mt:gt)??void 0,quad:Le,line:D,originalLine:U,curveBaseLine:nt,curvePoints:Ot,labelPoint:un,ridgePoints:N,outerSideMeans:b.outerSideMeans,outerSideQuads:b.outerSideQuads,distortionCorrected:e.distortionCorrected??!1})}}return(d=e.onProgress)==null||d.call(e,"Finalizing results...",.98),(g=e.onProgress)==null||g.call(e,"Finalizing results...",1),l}const Cc=n=>!n.blackLevels||n.blackLevels.length<4?null:[Number(n.blackLevels[0])||0,Number(n.blackLevels[1])||0,Number(n.blackLevels[2])||0,Number(n.blackLevels[3])||0],Gi=(n,t)=>{t instanceof ArrayBuffer&&(n.includes(t)||n.push(t))};self.onmessage=async n=>{var o,a;const{id:t,buffer:e,detect:i,options:r}=n.data,s=performance.now();try{const l=performance.now(),u=await io(e),h=performance.now()-l;let c=0,y=[];if(i&&!u.isXTrans){const p=u.isThreePlane?"three-plane":"rggb-raw",d=performance.now();y=Pc({width:u.width,height:u.height},u,{...r,sourceMode:p,forceRenderedMeasurement:!1,blackLevel:(r==null?void 0:r.blackLevel)??Cc(u),onProgress:(g,x)=>{self.postMessage({id:t,type:"progress",stage:g,progress:x})}}),c=performance.now()-d}const f=[];Gi(f,e),Gi(f,(o=u.data)==null?void 0:o.buffer),Gi(f,(a=u.floatData)==null?void 0:a.buffer),self.postMessage({id:t,type:"result",success:!0,raw:u,rawFileBuffer:e,measurements:y,timings:{decodeMs:h,detectMs:c,totalMs:performance.now()-s}},f)}catch(l){self.postMessage({id:t,type:"result",success:!1,error:(l==null?void 0:l.message)||String(l)})}};
