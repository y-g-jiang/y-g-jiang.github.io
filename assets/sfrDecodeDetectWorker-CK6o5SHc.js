var aa=Object.defineProperty;var oa=(n,t,e)=>t in n?aa(n,t,{enumerable:!0,configurable:!0,writable:!0,value:e}):n[t]=e;var yt=(n,t,e)=>oa(n,typeof t!="symbol"?t+"":t,e);function la(n){const[[t,e,i],[r,s,o],[a,l,u]]=n,h=t*(s*u-o*l)-e*(r*u-o*a)+i*(r*l-s*a);if(Math.abs(h)<1e-10)return null;const c=1/h;return[[(s*u-o*l)*c,(i*l-e*u)*c,(e*o-i*s)*c],[(o*a-r*u)*c,(t*u-i*a)*c,(i*r-t*o)*c],[(r*l-s*a)*c,(e*a-t*l)*c,(t*s-e*r)*c]]}function ca(n,t){const[[e,i,r],[s,o,a],[l,u,h]]=n,[c,m,d]=t;return[e*c+i*m+r*d,s*c+o*m+a*d,l*c+u*m+h*d]}class gn{constructor(t){yt(this,"size");yt(this,"isPowerOfTwo");yt(this,"_real");yt(this,"_imag");yt(this,"_scratch",null);yt(this,"_rev",null);yt(this,"_m",0);yt(this,"_internalFFT",null);yt(this,"_chirpReal",null);yt(this,"_chirpImag",null);yt(this,"_bReal",null);yt(this,"_bImag",null);yt(this,"_hanning",null);yt(this,"_windowSumSq",0);this.size=t,this.isPowerOfTwo=(t&t-1)===0&&t>0,this._real=new Float32Array(t),this._imag=new Float32Array(t),this.isPowerOfTwo?this.initRadix2():this.initBluestein()}initRadix2(){const t=this.size,e=Math.log2(t);this._rev=new Uint32Array(t);for(let i=0;i<t;i++){let r=0,s=i;for(let o=0;o<e;o++)r=r<<1|s&1,s>>>=1;this._rev[i]=r}}initBluestein(){const t=this.size;this._m=Math.pow(2,Math.ceil(Math.log2(2*t-1))),this._internalFFT=new gn(this._m),this._chirpReal=new Float32Array(t),this._chirpImag=new Float32Array(t);for(let r=0;r<t;r++){const s=-Math.PI*(r*r)/t;this._chirpReal[r]=Math.cos(s),this._chirpImag[r]=Math.sin(s)}const e=new Float32Array(this._m),i=new Float32Array(this._m);for(let r=0;r<t;r++)e[r]=this._chirpReal[r],i[r]=-this._chirpImag[r];for(let r=1;r<t;r++)e[this._m-r]=e[r],i[this._m-r]=i[r];this._internalFFT.transform(e,i),this._bReal=new Float32Array(this._internalFFT._real),this._bImag=new Float32Array(this._internalFFT._imag)}initHanning(){if(this._hanning)return;const t=this.size;this._hanning=new Float32Array(t);let e=0;for(let i=0;i<t;i++){const r=.5*(1-Math.cos(2*Math.PI*i/(t-1)));this._hanning[i]=r,e+=r*r}this._windowSumSq=e}transform(t,e){this.isPowerOfTwo?this.transformRadix2(t,e):this.transformBluestein(t,e)}transformRadix2(t,e){const i=this.size,r=this._rev,s=this._real,o=this._imag;if(t===s)for(let a=0;a<i;a++){const l=r[a];if(a<l){const u=s[a],h=o[a];s[a]=s[l],o[a]=o[l],s[l]=u,o[l]=h}}else for(let a=0;a<i;a++){const l=r[a];s[a]=t[l],o[a]=e?e[l]:0}for(let a=2;a<=i;a*=2){const l=a/2,u=-2*Math.PI/a,h=Math.cos(u),c=Math.sin(u);for(let m=0;m<i;m+=a){let d=1,p=0;for(let f=0;f<l;f++){const y=m+f,x=m+f+l,g=d*s[x]-p*o[x],b=d*o[x]+p*s[x],_=s[y],M=o[y];s[y]=_+g,o[y]=M+b,s[x]=_-g,o[x]=M-b;const w=d*h-p*c,P=d*c+p*h;d=w,p=P}}}}transformBluestein(t,e){const i=this.size,r=this._m,s=this._internalFFT,o=s._real,a=s._imag;o.fill(0),a.fill(0);for(let c=0;c<i;c++){const m=t[c],d=e?e[c]:0,p=this._chirpReal[c],f=this._chirpImag[c];o[c]=m*p-d*f,a[c]=m*f+d*p}s.transformRadix2(o,a);for(let c=0;c<r;c++){const m=s._real[c],d=s._imag[c],p=this._bReal[c],f=this._bImag[c];s._real[c]=m*p-d*f,s._imag[c]=m*f+d*p}const l=s._real,u=s._imag;for(let c=0;c<r;c++)u[c]=-u[c];s.transformRadix2(l,u);const h=1/r;for(let c=0;c<i;c++){const m=s._real[c]*h,d=-s._imag[c]*h,p=this._chirpReal[c],f=this._chirpImag[c];this._real[c]=m*p-d*f,this._imag[c]=m*f+d*p}}calculateSpectrum(t,e,i=!1){const r=this.size;let s=0;for(let h=0;h<r;h++)s+=t[h];const o=s/r;this._scratch||(this._scratch=new Float32Array(r));const a=this._scratch;if(i){this.initHanning();const h=this._hanning;for(let c=0;c<r;c++)a[c]=(t[c]-o)*h[c]}else for(let h=0;h<r;h++)a[h]=t[h]-o;this.transform(a);const l=e.length;let u=1/r;i&&this._windowSumSq>0&&(u=1/this._windowSumSq);for(let h=0;h<l;h++){const c=this._real[h],m=this._imag[h];e[h]+=(c*c+m*m)*u}}calculateSpectrumWindow(t,e,i,r=!1){const s=this.size;let o=0;for(let c=0;c<s;c++)o+=t[e+c];const a=o/s;this._scratch||(this._scratch=new Float32Array(s));const l=this._scratch;if(r){this.initHanning();const c=this._hanning;for(let m=0;m<s;m++)l[m]=(t[e+m]-a)*c[m]}else for(let c=0;c<s;c++)l[c]=t[e+c]-a;this.transform(l);const u=i.length;let h=1/s;r&&this._windowSumSq>0&&(h=1/this._windowSumSq);for(let c=0;c<u;c++){const m=this._real[c],d=this._imag[c];i[c]+=(m*m+d*d)*h}}}const ua={"Sony ILCE-7RM5":"0.82 -0.2976 -0.0719 -0.4296 1.2053 0.2532 -0.0429 0.1282 0.5774"};let Li=null;async function ha(n){return Li||(Li=(async()=>{if(typeof window.loadPyodide!="function")throw new Error("Pyodide missing: window.loadPyodide not found.");const t=await window.loadPyodide();return await t.loadPackage("numpy"),t})()),Li}let Ei=null,Tr=!1;async function tr(){var e;Ei||(Ei=(async()=>{const i=await import("./joraw2-Bb3_vNP4.js");if(typeof i.default!="function")throw new Error("JoRaw2 WASM import failed");const r=new URL("/assets/joraw2-3YkywkGx.wasm",import.meta.url).href;return i.default({locateFile(s,o){return s.endsWith("joraw2.wasm")?r:o+s}})})());const t=(await Ei).LibRaw;if(!t)throw new Error("JoRaw2 class not found");if(!Tr){const i=new t;try{if(typeof i.runtimeInfo!="function")throw new Error("JoRaw2 runtime identity is missing");const r=i.runtimeInfo();if((r==null?void 0:r.wrapper)!=="joraw2"||!String((r==null?void 0:r.librawVersion)||"").startsWith("0.22.2")||!(r!=null&&r.nikonHe)||!(r!=null&&r.nikonHeStar))throw new Error(`Unexpected JoRaw2 runtime: ${JSON.stringify(r)}`);Tr=!0}finally{(e=i.delete)==null||e.call(i)}}return t}var da=`#!/usr/bin/env python3
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
`,fa=`#!/usr/bin/env python3
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
`,pa=`#!/usr/bin/env python3
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
`;const Ir=512,ma=ua["Sony ILCE-7RM5"].split(/\s+/).map(Number).filter(Number.isFinite);let Ui=null;function ga(n){if(n.byteLength<8)return null;const t=n.getUint16(0,!1);return t===18761?!0:t===19789?!1:null}function ds(n,t,e){const i=Math.min(n.length,t+e);let r="";for(let s=t;s<i;s++){const o=n[s];if(o===0)break;r+=String.fromCharCode(o)}return r.trim()}function Gt(n,t,e,i,r){const s=e===1||e===2||e===7?1:e===3||e===8?2:e===4||e===9?4:0;if(!s)return[];const o=s*i,a=o<=4?r:n.getUint32(r,t);if(a<0||a+o>n.byteLength)return[];const l=[];for(let u=0;u<i;u++){const h=a+u*s;e===1||e===2||e===7?l.push(n.getUint8(h)):e===3?l.push(n.getUint16(h,t)):e===8?l.push(n.getInt16(h,t)):e===4?l.push(n.getUint32(h,t)):e===9&&l.push(n.getInt32(h,t))}return l}function Nr(n,t,e,i,r,s){if(i!==2||r<=0)return"";const o=r<=4?s:t.getUint32(s,e);return o<0||o>=n.length?"":ds(n,o,r)}function fs(n){const t=new Uint8Array(n),e=new DataView(n),i=ga(e);if(i===null||e.getUint16(2,i)!==42)return null;const s=c=>e.getUint16(c,i),o=c=>e.getUint32(c,i),a=[o(4)],l=new Set;let u="",h="";for(;a.length;){const c=a.pop();if(l.has(c)||c<=0||c+2>e.byteLength)continue;l.add(c);const m=s(c);if(c+2+m*12+4>e.byteLength)continue;const d=new Map;for(let _=0;_<m;_++){const M=c+2+_*12,w=s(M),P=s(M+2),C=o(M+4),v=M+8;d.set(w,{type:P,count:C,valueOffset:v})}const p=d.get(271),f=d.get(272);p&&!u&&(u=Nr(t,e,i,p.type,p.count,p.valueOffset)),f&&!h&&(h=Nr(t,e,i,f.type,f.count,f.valueOffset));const y=d.get(330);if(y){const _=Gt(e,i,y.type,y.count,y.valueOffset);for(const M of _)a.push(M)}const x=d.get(259),g=d.get(262);if(x&&g){const _=Gt(e,i,x.type,x.count,x.valueOffset)[0],M=Gt(e,i,g.type,g.count,g.valueOffset)[0];if(_===32766&&M===32803){const w=Gt(e,i,d.get(256).type,d.get(256).count,d.get(256).valueOffset)[0],P=Gt(e,i,d.get(257).type,d.get(257).count,d.get(257).valueOffset)[0],C=Gt(e,i,d.get(258).type,d.get(258).count,d.get(258).valueOffset)[0],v=Gt(e,i,d.get(273).type,d.get(273).count,d.get(273).valueOffset)[0],k=Gt(e,i,d.get(279).type,d.get(279).count,d.get(279).valueOffset)[0],A=d.get(33422)?Gt(e,i,d.get(33422).type,d.get(33422).count,d.get(33422).valueOffset):[0,1,1,2];d.get(29456)&&Gt(e,i,d.get(29456).type,d.get(29456).count,d.get(29456).valueOffset);const S=d.get(50717)?Gt(e,i,d.get(50717).type,d.get(50717).count,d.get(50717).valueOffset)[0]:16383,F=d.get(50719)?Gt(e,i,d.get(50719).type,d.get(50719).count,d.get(50719).valueOffset):[],T=d.get(50720)?Gt(e,i,d.get(50720).type,d.get(50720).count,d.get(50720).valueOffset):[];if(v+Ir+16>t.length||v+k>t.length)return null;const R=v+Ir,E=ds(t,R,4),O=t[R+8]<<8|t[R+9],U=t[R+10]<<8|t[R+11],I=t[R+12]<<8|t[R+13],D=t[R+14]<<8|t[R+15],B=I>>4&63,V=D>>13,Q=D>>10&3,z=U*2,Y=k>=4?(t[v]|t[v+1]<<8|t[v+2]<<16|t[v+3]<<24)>>>0:0,W=O===w&&z===P;let K=!1;if(Y>=1&&Y<=16&&k>=8+Y*24){const Z=new Map,L=new Map;let G=!0;for(let X=0;X<Y;X++){const j=v+8+X*24,$=e.getUint32(j+8,!0),et=e.getUint32(j+12,!0),nt=e.getUint32(j+16,!0),ct=e.getUint32(j+20,!0);if(!nt||!ct||$+nt>w||et+ct>P){G=!1;break}const at=Z.get(et);if(at!==void 0&&at!==ct){G=!1;break}Z.set(et,ct),L.set(et,(L.get(et)||0)+nt)}if(G){const X=Array.from(Z.keys()).sort(($,et)=>$-et);let j=0;for(const $ of X){if($!==j||L.get($)!==w){G=!1;break}j+=Z.get($)}K=G&&j===P}}const it=Y>=1&&Y<=16&&O>0&&z>0&&w%O===0&&P%z===0&&Y===w/O*(P/z);if(E!=="A000"&&E!=="0000"||!W&&!it&&!K||B!==16||V!==3||Q!==3)return null;const J=[1024,1024,1024,1024],ot=A.length>=4?A.slice(0,4).map(Z=>Z===0?"R":Z===2?"B":Z===1?"G":"?").join(""):"";return{width:w,height:P,bitsPerSample:C,compression:_,photometric:M,blackLevel:J,whiteLevel:Number(S||16383),cfaPattern:ot,defaultCropOrigin:F.length>=2?[Number(F[0]),Number(F[1])]:void 0,defaultCropSize:T.length>=2?[Number(T[0]),Number(T[1])]:void 0,make:u||"SONY",model:h||"ILCE-7M5"}}}const b=o(c+2+m*12);b&&a.push(b)}return null}async function ya(n){return Ui||(Ui=(async()=>{const t=await ha();return t.__jtrSonyCrawHqDecoderReady||(await t.FS.mkdirTree("/sony_craw_hq"),await t.FS.writeFile("/sony_craw_hq/llvc3_bitstream_probe.py",da),await t.FS.writeFile("/sony_craw_hq/llvc3_entropy.py",fa),await t.FS.writeFile("/sony_craw_hq/llvc3_math.py",pa),await t.runPythonAsync(`
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
`),t.__jtrSonyCrawHqDecoderReady=!0),t})()),Ui}function xa(n){return fs(n)}function ba(n){return n==="ILCE-7M5"?ma:null}function ps(n,t,e,i){var l;if(n.length!==t.width*t.height)throw new Error(`Sony cRAW HQ decoded size mismatch: got ${n.length}, expected ${t.width*t.height}`);const r=t.model||"ILCE-7M5",s=r.startsWith("Sony ")?r:`Sony ${r}`,o=ba(r),a={...i||{},make:t.make||(i==null?void 0:i.camera_make)||"SONY",model:r,camera_make:t.make||(i==null?void 0:i.camera_make)||"SONY",camera_model:r,UniqueCameraModel:s,sourceFormat:e==="joraw2-wasm"?"Sony cRAW HQ / LLVC3 (JoRaw2 WASM)":"Sony cRAW HQ / LLVC3",sonyCrawHq:{...t,decodeBackend:e},color_desc:t.cfaPattern,black_level_per_channel:t.blackLevel,white_level:t.whiteLevel,color_matrix:o&&o.length===9?o:void 0,idata:{filters:2492765332,colors:3},color_data:{...(i==null?void 0:i.color_data)||{},black:1024,cblack_rawpy_style:t.blackLevel,dng_levels:{...((l=i==null?void 0:i.color_data)==null?void 0:l.dng_levels)||{},dng_cblack:t.blackLevel,dng_whitelevel:t.whiteLevel}}};return{data:n,width:t.width,height:t.height,bayerPattern:t.cfaPattern,blackLevels:t.blackLevel,whiteLevel:t.whiteLevel,metadata:a,isThreePlane:!1,isXTrans:!1}}async function _a(n,t,e){const i=typeof performance<"u"?performance.now():Date.now(),r=await tr(),s=typeof performance<"u"?performance.now():Date.now(),o=new r;try{const a=new Uint8Array(n);await o.open(a,{});const l=typeof performance<"u"?performance.now():Date.now();let u=null;try{u=await o.metadata(!0)}catch(y){console.warn("[Sony cRAW HQ] fast WASM metadata read failed",y)}const h=typeof performance<"u"?performance.now():Date.now(),c=o.getRawImage(),m=typeof performance<"u"?performance.now():Date.now();if(!c||!c.data)throw new Error("Sony cRAW HQ LibRaw WASM returned no raw image");const d=c.data instanceof Uint16Array?c.data:new Uint16Array(c.data.buffer,c.data.byteOffset||0,c.data.byteLength/2),p=typeof performance<"u"?performance.now():Date.now();if(c.width!==t.width||c.height!==t.height)throw new Error(`Sony cRAW HQ LibRaw WASM dimensions mismatch: got ${c.width}x${c.height}, expected ${t.width}x${t.height}`);const f=typeof performance<"u"?performance.now():Date.now();return console.info("[Sony cRAW HQ] fast decode timings",{width:t.width,height:t.height,backend:"joraw2-wasm",wasmReadyMs:Math.round(s-i),openMs:Math.round(l-s),metadataMs:Math.round(h-l),unpackMs:Math.round(m-h),copyMs:Math.round(p-m),totalMs:Math.round(f-i)}),{rawImageData:ps(d,t,"joraw2-wasm",u),info:t}}finally{typeof o.delete=="function"?o.delete():typeof o.close=="function"&&o.close()}}async function wa(n,t,e){const i=typeof performance<"u"?performance.now():Date.now(),r=await ya(),s=typeof performance<"u"?performance.now():Date.now(),o=new Uint8Array(n),a=await fetch(new URL("/assets/sony_llvc3_static_lut4096_padded_u16-FsVBk-IV.bin",import.meta.url));if(!a.ok)throw new Error(`Failed to load Sony LLVC3 sample LUT: HTTP ${a.status}`);const l=new Uint8Array(await a.arrayBuffer()),u=typeof performance<"u"?performance.now():Date.now();r.globals.set("jtr_sony_arw_bytes",o),r.globals.set("jtr_sony_lut_bytes",l);const h=await r.runPythonAsync("jtr_decode_sony_craw_hq(jtr_sony_arw_bytes.to_py(), jtr_sony_lut_bytes.to_py())"),c=typeof performance<"u"?performance.now():Date.now(),m=h.toJs();typeof h.destroy=="function"&&h.destroy(),r.globals.delete("jtr_sony_arw_bytes"),r.globals.delete("jtr_sony_lut_bytes");const d=new Uint8Array(m.byteLength);d.set(m);const p=new Uint16Array(d.buffer),f=typeof performance<"u"?performance.now():Date.now(),y=typeof performance<"u"?performance.now():Date.now();return console.info("[Sony cRAW HQ] decode timings",{width:t.width,height:t.height,backend:"pyodide",pyodideReadyMs:Math.round(s-i),lutLoadMs:Math.round(u-s),llvc3DecodeMs:Math.round(c-u),copyMs:Math.round(f-c),totalMs:Math.round(y-i)}),{rawImageData:ps(p,t,"pyodide"),info:t}}async function Ma(n,t){const e=fs(n);if(!e)return null;try{return await _a(n,e,t)}catch(i){return console.warn("[Sony cRAW HQ] fast WASM decode failed; falling back to Pyodide",i),wa(n,e)}}async function Sa(n,t){return Ma(n,t)}const Pa=["RGGB","BGGR","GRBG","GBRG"],va=new Set(Pa);function Ca(n){if(!n||typeof n!="object"||Array.isArray(n)||ArrayBuffer.isView(n))return n;const t=n;return t.value??t.values??t.description??n}function $e(n){const t=Ca(n);let e="";if(typeof t=="string")e=t.toUpperCase().replace(/[^RGB012]/g,"");else if(typeof t=="number")e=String(t);else if(Array.isArray(t)||ArrayBuffer.isView(t))e=Array.from(t).map(String).join("");else return null;return/^[012]{4}$/.test(e)&&(e=e.replace(/0/g,"R").replace(/1/g,"G").replace(/2/g,"B")),va.has(e)?e:null}function Rr(n){var t;return n?[n.cfa_pattern,n.cfaPattern,n.BayerPattern,n.CFAPattern2,n.CFAPattern,(t=n.idata)==null?void 0:t.cfa_pattern]:[]}function Fa(n,t,e){const i=$e(n.bayerPattern),r=n.bayerPatternSource==="manual";if(i&&!r)return{pattern:i,source:n.bayerPatternSource||"decoder"};if(!r)for(const o of Rr(n.metadata)){const a=$e(o);if(a)return{pattern:a,source:"libraw-metadata"}}for(const o of Rr(t)){const a=$e(o);if(a)return{pattern:a,source:"metadata"}}const s=$e(e)||(r?i:null);return s?{pattern:s,source:"manual"}:{pattern:null,source:null}}function ka(n,t){return n.bayerPattern=t.pattern||"",n.bayerPatternSource=t.source||void 0,t.pattern&&t.source!=="manual"&&(n.metadata={...n.metadata||{},cfa_pattern:t.pattern,cfaPattern:t.pattern}),t}function ms(n,t){return(t&1)<<1|n&1}function Aa(n,t,e){const i=$e(n);if(!i)throw new Error("Bayer CFA pattern is unresolved.");return i[ms(t,e)]}const ce=n=>{const t=Number(n);return Number.isFinite(t)?Math.max(0,t):0};function Ke(n){if(!n||typeof n.length!="number")return null;const t=Array.from(n);return t.length<4?null:[ce(t[0]),ce(t[1]),ce(t[2]),ce(t[3])]}function Lr(n){if(!n||n.source!=="libraw")return null;const t=Ke(n.channelOffsets),e=Ke(n.channelLevels),i=Ke(n.siteColorIndices),r=Ke(n.siteBaseLevels),s=Ke(n.siteLevels);if(!t||!e||!i||!r||!s)return null;const o=Math.max(0,Math.floor(ce(n.repeatRows))),a=Math.max(0,Math.floor(ce(n.repeatCols))),l=o*a,u=n.repeatValues&&typeof n.repeatValues.length=="number"?Array.from(n.repeatValues,ce):[],h=l>0&&l<=4098&&u.length>=l;return{source:"libraw",common:ce(n.common),channelOffsets:t,channelLevels:e,siteColorIndices:i,siteBaseLevels:r,siteLevels:s,repeatRows:h?o:0,repeatCols:h?a:0,repeatOriginY:Math.max(0,Math.floor(ce(n.repeatOriginY))),repeatOriginX:Math.max(0,Math.floor(ce(n.repeatOriginX))),repeatValues:h?u.slice(0,l):[]}}var fi=typeof self<"u"?self:global;const Rn=typeof navigator<"u",Ta=Rn&&typeof HTMLImageElement>"u",ei=!(typeof global>"u"||typeof process>"u"||!process.versions||!process.versions.node),pi=fi.Buffer,Hn=fi.BigInt,mi=!!pi,Ia=n=>n;function ni(n,t=Ia){if(ei)try{return typeof require=="function"?Promise.resolve(t(require(n))):import(n).then(t)}catch{console.warn(`Couldn't load ${n}`)}}let er=fi.fetch;const Na=n=>er=n;if(!fi.fetch){const n=ni("http",i=>i),t=ni("https",i=>i),e=(i,{headers:r}={})=>new Promise(async(s,o)=>{let{port:a,hostname:l,pathname:u,protocol:h,search:c}=new URL(i);const m={method:"GET",hostname:l,path:encodeURI(u)+c,headers:r};a!==""&&(m.port=Number(a));const d=(h==="https:"?await t:await n).request(m,p=>{if(p.statusCode===301||p.statusCode===302){let f=new URL(p.headers.location,i).toString();return e(f,{headers:r}).then(s).catch(o)}s({status:p.statusCode,arrayBuffer:()=>new Promise(f=>{let y=[];p.on("data",x=>y.push(x)),p.on("end",()=>f(Buffer.concat(y)))})})});d.on("error",o),d.end()});Na(e)}function tt(n,t,e){return t in n?Object.defineProperty(n,t,{value:e,enumerable:!0,configurable:!0,writable:!0}):n[t]=e,n}const ii=n=>gs(n)?void 0:n,Ra=n=>n!==void 0;function gs(n){return n===void 0||(n instanceof Map?n.size===0:Object.values(n).filter(Ra).length===0)}function Pt(n){let t=new Error(n);throw delete t.stack,t}function Je(n){return(n=function(t){for(;t.endsWith("\0");)t=t.slice(0,-1);return t}(n).trim())===""?void 0:n}function Hi(n){let t=function(e){let i=0;return e.ifd0.enabled&&(i+=1024),e.exif.enabled&&(i+=2048),e.makerNote&&(i+=2048),e.userComment&&(i+=1024),e.gps.enabled&&(i+=512),e.interop.enabled&&(i+=100),e.ifd1.enabled&&(i+=1024),i+2048}(n);return n.jfif.enabled&&(t+=50),n.xmp.enabled&&(t+=2e4),n.iptc.enabled&&(t+=14e3),n.icc.enabled&&(t+=6e3),t}const ji=n=>String.fromCharCode.apply(null,n),Er=typeof TextDecoder<"u"?new TextDecoder("utf-8"):void 0;function ys(n){return Er?Er.decode(n):mi?Buffer.from(n).toString("utf8"):decodeURIComponent(escape(ji(n)))}class zt{static from(t,e){return t instanceof this&&t.le===e?t:new zt(t,void 0,void 0,e)}constructor(t,e=0,i,r){if(typeof r=="boolean"&&(this.le=r),Array.isArray(t)&&(t=new Uint8Array(t)),t===0)this.byteOffset=0,this.byteLength=0;else if(t instanceof ArrayBuffer){i===void 0&&(i=t.byteLength-e);let s=new DataView(t,e,i);this._swapDataView(s)}else if(t instanceof Uint8Array||t instanceof DataView||t instanceof zt){i===void 0&&(i=t.byteLength-e),(e+=t.byteOffset)+i>t.byteOffset+t.byteLength&&Pt("Creating view outside of available memory in ArrayBuffer");let s=new DataView(t.buffer,e,i);this._swapDataView(s)}else if(typeof t=="number"){let s=new DataView(new ArrayBuffer(t));this._swapDataView(s)}else Pt("Invalid input argument for BufferView: "+t)}_swapArrayBuffer(t){this._swapDataView(new DataView(t))}_swapBuffer(t){this._swapDataView(new DataView(t.buffer,t.byteOffset,t.byteLength))}_swapDataView(t){this.dataView=t,this.buffer=t.buffer,this.byteOffset=t.byteOffset,this.byteLength=t.byteLength}_lengthToEnd(t){return this.byteLength-t}set(t,e,i=zt){return t instanceof DataView||t instanceof zt?t=new Uint8Array(t.buffer,t.byteOffset,t.byteLength):t instanceof ArrayBuffer&&(t=new Uint8Array(t)),t instanceof Uint8Array||Pt("BufferView.set(): Invalid data argument."),this.toUint8().set(t,e),new i(this,e,t.byteLength)}subarray(t,e){return e=e||this._lengthToEnd(t),new zt(this,t,e)}toUint8(){return new Uint8Array(this.buffer,this.byteOffset,this.byteLength)}getUint8Array(t,e){return new Uint8Array(this.buffer,this.byteOffset+t,e)}getString(t=0,e=this.byteLength){return ys(this.getUint8Array(t,e))}getLatin1String(t=0,e=this.byteLength){let i=this.getUint8Array(t,e);return ji(i)}getUnicodeString(t=0,e=this.byteLength){const i=[];for(let r=0;r<e&&t+r<this.byteLength;r+=2)i.push(this.getUint16(t+r));return ji(i)}getInt8(t){return this.dataView.getInt8(t)}getUint8(t){return this.dataView.getUint8(t)}getInt16(t,e=this.le){return this.dataView.getInt16(t,e)}getInt32(t,e=this.le){return this.dataView.getInt32(t,e)}getUint16(t,e=this.le){return this.dataView.getUint16(t,e)}getUint32(t,e=this.le){return this.dataView.getUint32(t,e)}getFloat32(t,e=this.le){return this.dataView.getFloat32(t,e)}getFloat64(t,e=this.le){return this.dataView.getFloat64(t,e)}getFloat(t,e=this.le){return this.dataView.getFloat32(t,e)}getDouble(t,e=this.le){return this.dataView.getFloat64(t,e)}getUintBytes(t,e,i){switch(e){case 1:return this.getUint8(t,i);case 2:return this.getUint16(t,i);case 4:return this.getUint32(t,i);case 8:return this.getUint64&&this.getUint64(t,i)}}getUint(t,e,i){switch(e){case 8:return this.getUint8(t,i);case 16:return this.getUint16(t,i);case 32:return this.getUint32(t,i);case 64:return this.getUint64&&this.getUint64(t,i)}}toString(t){return this.dataView.toString(t,this.constructor.name)}ensureChunk(){}}function qi(n,t){Pt(`${n} '${t}' was not loaded, try using full build of exifr.`)}class nr extends Map{constructor(t){super(),this.kind=t}get(t,e){return this.has(t)||qi(this.kind,t),e&&(t in e||function(i,r){Pt(`Unknown ${i} '${r}'.`)}(this.kind,t),e[t].enabled||qi(this.kind,t)),super.get(t)}keyList(){return Array.from(this.keys())}}var he=new nr("file parser"),Mt=new nr("segment parser"),fe=new nr("file reader");function La(n,t){return typeof n=="string"?Ur(n,t):Rn&&!Ta&&n instanceof HTMLImageElement?Ur(n.src,t):n instanceof Uint8Array||n instanceof ArrayBuffer||n instanceof DataView?new zt(n):Rn&&n instanceof Blob?Qi(n,t,"blob",an):void Pt("Invalid input argument")}function Ur(n,t){return(e=n).startsWith("data:")||e.length>1e4?Ki(n,t,"base64"):ei&&n.includes("://")?Qi(n,t,"url",sn):ei?Ki(n,t,"fs"):Rn?Qi(n,t,"url",sn):void Pt("Invalid input argument");var e}async function Qi(n,t,e,i){return fe.has(e)?Ki(n,t,e):i?async function(r,s){let o=await s(r);return new zt(o)}(n,i):void Pt(`Parser ${e} is not loaded`)}async function Ki(n,t,e){let i=new(fe.get(e))(n,t);return await i.read(),i}const sn=n=>er(n).then(t=>t.arrayBuffer()),an=n=>new Promise((t,e)=>{let i=new FileReader;i.onloadend=()=>t(i.result||new ArrayBuffer),i.onerror=e,i.readAsArrayBuffer(n)});class Ea extends Map{get tagKeys(){return this.allKeys||(this.allKeys=Array.from(this.keys())),this.allKeys}get tagValues(){return this.allValues||(this.allValues=Array.from(this.values())),this.allValues}}function wt(n,t,e){let i=new Ea;for(let[r,s]of e)i.set(r,s);if(Array.isArray(t))for(let r of t)n.set(r,i);else n.set(t,i);return i}function on(n,t,e){let i,r=n.get(t);for(i of e)r.set(i[0],i[1])}const kt=new Map,$t=new Map,Ne=new Map,Ce=["chunked","firstChunkSize","firstChunkSizeNode","firstChunkSizeBrowser","chunkSize","chunkLimit"],yn=["jfif","xmp","icc","iptc","ihdr"],ln=["tiff",...yn],bt=["ifd0","ifd1","exif","gps","interop"],Fe=[...ln,...bt],ke=["makerNote","userComment"],xn=["translateKeys","translateValues","reviveValues","multiSegment"],Ae=[...xn,"sanitize","mergeOutput","silentErrors"];class xs{get translate(){return this.translateKeys||this.translateValues||this.reviveValues}}class Fn extends xs{get needed(){return this.enabled||this.deps.size>0}constructor(t,e,i,r){if(super(),tt(this,"enabled",!1),tt(this,"skip",new Set),tt(this,"pick",new Set),tt(this,"deps",new Set),tt(this,"translateKeys",!1),tt(this,"translateValues",!1),tt(this,"reviveValues",!1),this.key=t,this.enabled=e,this.parse=this.enabled,this.applyInheritables(r),this.canBeFiltered=bt.includes(t),this.canBeFiltered&&(this.dict=kt.get(t)),i!==void 0)if(Array.isArray(i))this.parse=this.enabled=!0,this.canBeFiltered&&i.length>0&&this.translateTagSet(i,this.pick);else if(typeof i=="object"){if(this.enabled=!0,this.parse=i.parse!==!1,this.canBeFiltered){let{pick:s,skip:o}=i;s&&s.length>0&&this.translateTagSet(s,this.pick),o&&o.length>0&&this.translateTagSet(o,this.skip)}this.applyInheritables(i)}else i===!0||i===!1?this.parse=this.enabled=i:Pt(`Invalid options argument: ${i}`)}applyInheritables(t){let e,i;for(e of xn)i=t[e],i!==void 0&&(this[e]=i)}translateTagSet(t,e){if(this.dict){let i,r,{tagKeys:s,tagValues:o}=this.dict;for(i of t)typeof i=="string"?(r=o.indexOf(i),r===-1&&(r=s.indexOf(Number(i))),r!==-1&&e.add(Number(s[r]))):e.add(i)}else for(let i of t)e.add(i)}finalizeFilters(){!this.enabled&&this.deps.size>0?(this.enabled=!0,ri(this.pick,this.deps)):this.enabled&&this.pick.size>0&&ri(this.pick,this.deps)}}var Lt={jfif:!1,tiff:!0,xmp:!1,icc:!1,iptc:!1,ifd0:!0,ifd1:!1,exif:!0,gps:!0,interop:!1,ihdr:void 0,makerNote:!1,userComment:!1,multiSegment:!1,skip:[],pick:[],translateKeys:!0,translateValues:!0,reviveValues:!0,sanitize:!0,mergeOutput:!0,silentErrors:!0,chunked:!0,firstChunkSize:void 0,firstChunkSizeNode:512,firstChunkSizeBrowser:65536,chunkSize:65536,chunkLimit:5},Dr=new Map;class cn extends xs{static useCached(t){let e=Dr.get(t);return e!==void 0||(e=new this(t),Dr.set(t,e)),e}constructor(t){super(),t===!0?this.setupFromTrue():t===void 0?this.setupFromUndefined():Array.isArray(t)?this.setupFromArray(t):typeof t=="object"?this.setupFromObject(t):Pt(`Invalid options argument ${t}`),this.firstChunkSize===void 0&&(this.firstChunkSize=Rn?this.firstChunkSizeBrowser:this.firstChunkSizeNode),this.mergeOutput&&(this.ifd1.enabled=!1),this.filterNestedSegmentTags(),this.traverseTiffDependencyTree(),this.checkLoadedPlugins()}setupFromUndefined(){let t;for(t of Ce)this[t]=Lt[t];for(t of Ae)this[t]=Lt[t];for(t of ke)this[t]=Lt[t];for(t of Fe)this[t]=new Fn(t,Lt[t],void 0,this)}setupFromTrue(){let t;for(t of Ce)this[t]=Lt[t];for(t of Ae)this[t]=Lt[t];for(t of ke)this[t]=!0;for(t of Fe)this[t]=new Fn(t,!0,void 0,this)}setupFromArray(t){let e;for(e of Ce)this[e]=Lt[e];for(e of Ae)this[e]=Lt[e];for(e of ke)this[e]=Lt[e];for(e of Fe)this[e]=new Fn(e,!1,void 0,this);this.setupGlobalFilters(t,void 0,bt)}setupFromObject(t){let e;for(e of(bt.ifd0=bt.ifd0||bt.image,bt.ifd1=bt.ifd1||bt.thumbnail,Object.assign(this,t),Ce))this[e]=Di(t[e],Lt[e]);for(e of Ae)this[e]=Di(t[e],Lt[e]);for(e of ke)this[e]=Di(t[e],Lt[e]);for(e of ln)this[e]=new Fn(e,Lt[e],t[e],this);for(e of bt)this[e]=new Fn(e,Lt[e],t[e],this.tiff);this.setupGlobalFilters(t.pick,t.skip,bt,Fe),t.tiff===!0?this.batchEnableWithBool(bt,!0):t.tiff===!1?this.batchEnableWithUserValue(bt,t):Array.isArray(t.tiff)?this.setupGlobalFilters(t.tiff,void 0,bt):typeof t.tiff=="object"&&this.setupGlobalFilters(t.tiff.pick,t.tiff.skip,bt)}batchEnableWithBool(t,e){for(let i of t)this[i].enabled=e}batchEnableWithUserValue(t,e){for(let i of t){let r=e[i];this[i].enabled=r!==!1&&r!==void 0}}setupGlobalFilters(t,e,i,r=i){if(t&&t.length){for(let o of r)this[o].enabled=!1;let s=Br(t,i);for(let[o,a]of s)ri(this[o].pick,a),this[o].enabled=!0}else if(e&&e.length){let s=Br(e,i);for(let[o,a]of s)ri(this[o].skip,a)}}filterNestedSegmentTags(){let{ifd0:t,exif:e,xmp:i,iptc:r,icc:s}=this;this.makerNote?e.deps.add(37500):e.skip.add(37500),this.userComment?e.deps.add(37510):e.skip.add(37510),i.enabled||t.skip.add(700),r.enabled||t.skip.add(33723),s.enabled||t.skip.add(34675)}traverseTiffDependencyTree(){let{ifd0:t,exif:e,gps:i,interop:r}=this;r.needed&&(e.deps.add(40965),t.deps.add(40965)),e.needed&&t.deps.add(34665),i.needed&&t.deps.add(34853),this.tiff.enabled=bt.some(s=>this[s].enabled===!0)||this.makerNote||this.userComment;for(let s of bt)this[s].finalizeFilters()}get onlyTiff(){return!yn.map(t=>this[t].enabled).some(t=>t===!0)&&this.tiff.enabled}checkLoadedPlugins(){for(let t of ln)this[t].enabled&&!Mt.has(t)&&qi("segment parser",t)}}function Br(n,t){let e,i,r,s,o=[];for(r of t){for(s of(e=kt.get(r),i=[],e))(n.includes(s[0])||n.includes(s[1]))&&i.push(s[0]);i.length&&o.push([r,i])}return o}function Di(n,t){return n!==void 0?n:t!==void 0?t:void 0}function ri(n,t){for(let e of t)n.add(e)}tt(cn,"default",Lt);class Re{constructor(t){tt(this,"parsers",{}),tt(this,"output",{}),tt(this,"errors",[]),tt(this,"pushToErrors",e=>this.errors.push(e)),this.options=cn.useCached(t)}async read(t){this.file=await La(t,this.options)}setup(){if(this.fileParser)return;let{file:t}=this,e=t.getUint16(0);for(let[i,r]of he)if(r.canHandle(t,e))return this.fileParser=new r(this.options,this.file,this.parsers),t[i]=!0;this.file.close&&this.file.close(),Pt("Unknown file format")}async parse(){let{output:t,errors:e}=this;return this.setup(),this.options.silentErrors?(await this.executeParsers().catch(this.pushToErrors),e.push(...this.fileParser.errors)):await this.executeParsers(),this.file.close&&this.file.close(),this.options.silentErrors&&e.length>0&&(t.errors=e),ii(t)}async executeParsers(){let{output:t}=this;await this.fileParser.parse();let e=Object.values(this.parsers).map(async i=>{let r=await i.parse();i.assignToOutput(t,r)});this.options.silentErrors&&(e=e.map(i=>i.catch(this.pushToErrors))),await Promise.all(e)}async extractThumbnail(){this.setup();let{options:t,file:e}=this,i=Mt.get("tiff",t);var r;if(e.tiff?r={start:0,type:"tiff"}:e.jpeg&&(r=await this.fileParser.getOrFindSegment("tiff")),r===void 0)return;let s=await this.fileParser.ensureSegmentChunk(r),o=this.parsers.tiff=new i(s,t,e),a=await o.extractThumbnail();return e.close&&e.close(),a}}async function gi(n,t){let e=new Re(t);return await e.read(n),e.parse()}var Ua=Object.freeze({__proto__:null,parse:gi,Exifr:Re,fileParsers:he,segmentParsers:Mt,fileReaders:fe,tagKeys:kt,tagValues:$t,tagRevivers:Ne,createDictionary:wt,extendDictionary:on,fetchUrlAsArrayBuffer:sn,readBlobAsArrayBuffer:an,chunkedProps:Ce,otherSegments:yn,segments:ln,tiffBlocks:bt,segmentsAndBlocks:Fe,tiffExtractables:ke,inheritables:xn,allFormatters:Ae,Options:cn});class yi{constructor(t,e,i){tt(this,"errors",[]),tt(this,"ensureSegmentChunk",async r=>{let s=r.start,o=r.size||65536;if(this.file.chunked)if(this.file.available(s,o))r.chunk=this.file.subarray(s,o);else try{r.chunk=await this.file.readChunk(s,o)}catch(a){Pt(`Couldn't read segment: ${JSON.stringify(r)}. ${a.message}`)}else this.file.byteLength>s+o?r.chunk=this.file.subarray(s,o):r.size===void 0?r.chunk=this.file.subarray(s):Pt("Segment unreachable: "+JSON.stringify(r));return r.chunk}),this.extendOptions&&this.extendOptions(t),this.options=t,this.file=e,this.parsers=i}injectSegment(t,e){this.options[t].enabled&&this.createParser(t,e)}createParser(t,e){let i=new(Mt.get(t))(e,this.options,this.file);return this.parsers[t]=i}createParsers(t){for(let e of t){let{type:i,chunk:r}=e,s=this.options[i];if(s&&s.enabled){let o=this.parsers[i];o&&o.append||o||this.createParser(i,r)}}}async readSegments(t){let e=t.map(this.ensureSegmentChunk);await Promise.all(e)}}class Kt{static findPosition(t,e){let i=t.getUint16(e+2)+2,r=typeof this.headerLength=="function"?this.headerLength(t,e,i):this.headerLength,s=e+r,o=i-r;return{offset:e,length:i,headerLength:r,start:s,size:o,end:s+o}}static parse(t,e={}){return new this(t,new cn({[this.type]:e}),t).parse()}normalizeInput(t){return t instanceof zt?t:new zt(t)}constructor(t,e={},i){tt(this,"errors",[]),tt(this,"raw",new Map),tt(this,"handleError",r=>{if(!this.options.silentErrors)throw r;this.errors.push(r.message)}),this.chunk=this.normalizeInput(t),this.file=i,this.type=this.constructor.type,this.globalOptions=this.options=e,this.localOptions=e[this.type],this.canTranslate=this.localOptions&&this.localOptions.translate}translate(){this.canTranslate&&(this.translated=this.translateBlock(this.raw,this.type))}get output(){return this.translated?this.translated:this.raw?Object.fromEntries(this.raw):void 0}translateBlock(t,e){let i=Ne.get(e),r=$t.get(e),s=kt.get(e),o=this.options[e],a=o.reviveValues&&!!i,l=o.translateValues&&!!r,u=o.translateKeys&&!!s,h={};for(let[c,m]of t)a&&i.has(c)?m=i.get(c)(m):l&&r.has(c)&&(m=this.translateValue(m,r.get(c))),u&&s.has(c)&&(c=s.get(c)||c),h[c]=m;return h}translateValue(t,e){return e[t]||e.DEFAULT||t}assignToOutput(t,e){this.assignObjectToOutput(t,this.constructor.type,e)}assignObjectToOutput(t,e,i){if(this.globalOptions.mergeOutput)return Object.assign(t,i);t[e]?Object.assign(t[e],i):t[e]=i}}tt(Kt,"headerLength",4),tt(Kt,"type",void 0),tt(Kt,"multiSegment",!1),tt(Kt,"canHandle",()=>!1);function Da(n){return n===192||n===194||n===196||n===219||n===221||n===218||n===254}function Ba(n){return n>=224&&n<=239}function Oa(n,t,e){for(let[i,r]of Mt)if(r.canHandle(n,t,e))return i}class Or extends yi{constructor(...t){super(...t),tt(this,"appSegments",[]),tt(this,"jpegSegments",[]),tt(this,"unknownSegments",[])}static canHandle(t,e){return e===65496}async parse(){await this.findAppSegments(),await this.readSegments(this.appSegments),this.mergeMultiSegments(),this.createParsers(this.mergedAppSegments||this.appSegments)}setupSegmentFinderArgs(t){t===!0?(this.findAll=!0,this.wanted=new Set(Mt.keyList())):(t=t===void 0?Mt.keyList().filter(e=>this.options[e].enabled):t.filter(e=>this.options[e].enabled&&Mt.has(e)),this.findAll=!1,this.remaining=new Set(t),this.wanted=new Set(t)),this.unfinishedMultiSegment=!1}async findAppSegments(t=0,e){this.setupSegmentFinderArgs(e);let{file:i,findAll:r,wanted:s,remaining:o}=this;if(!r&&this.file.chunked&&(r=Array.from(s).some(a=>{let l=Mt.get(a),u=this.options[a];return l.multiSegment&&u.multiSegment}),r&&await this.file.readWhole()),t=this.findAppSegmentsInRange(t,i.byteLength),!this.options.onlyTiff&&i.chunked){let a=!1;for(;o.size>0&&!a&&(i.canReadNextChunk||this.unfinishedMultiSegment);){let{nextChunkOffset:l}=i,u=this.appSegments.some(h=>!this.file.available(h.offset||h.start,h.length||h.size));if(a=t>l&&!u?!await i.readNextChunk(t):!await i.readNextChunk(l),(t=this.findAppSegmentsInRange(t,i.byteLength))===void 0)return}}}findAppSegmentsInRange(t,e){e-=2;let i,r,s,o,a,l,{file:u,findAll:h,wanted:c,remaining:m,options:d}=this;for(;t<e;t++)if(u.getUint8(t)===255){if(i=u.getUint8(t+1),Ba(i)){if(r=u.getUint16(t+2),s=Oa(u,t,r),s&&c.has(s)&&(o=Mt.get(s),a=o.findPosition(u,t),l=d[s],a.type=s,this.appSegments.push(a),!h&&(o.multiSegment&&l.multiSegment?(this.unfinishedMultiSegment=a.chunkNumber<a.chunkCount,this.unfinishedMultiSegment||m.delete(s)):m.delete(s),m.size===0)))break;d.recordUnknownSegments&&(a=Kt.findPosition(u,t),a.marker=i,this.unknownSegments.push(a)),t+=r+1}else if(Da(i)){if(r=u.getUint16(t+2),i===218&&d.stopAfterSos!==!1)return;d.recordJpegSegments&&this.jpegSegments.push({offset:t,length:r,marker:i}),t+=r+1}}return t}mergeMultiSegments(){if(!this.appSegments.some(e=>e.multiSegment))return;let t=function(e,i){let r,s,o,a=new Map;for(let l=0;l<e.length;l++)r=e[l],s=r[i],a.has(s)?o=a.get(s):a.set(s,o=[]),o.push(r);return Array.from(a)}(this.appSegments,"type");this.mergedAppSegments=t.map(([e,i])=>{let r=Mt.get(e,this.options);return r.handleMultiSegments?{type:e,chunk:r.handleMultiSegments(i)}:i[0]})}getSegment(t){return this.appSegments.find(e=>e.type===t)}async getOrFindSegment(t){let e=this.getSegment(t);return e===void 0&&(await this.findAppSegments(0,[t]),e=this.getSegment(t)),e}}tt(Or,"type","jpeg"),he.set("jpeg",Or);const za=[void 0,1,1,2,4,8,1,1,2,4,8,4,8,4];class Va extends Kt{parseHeader(){var t=this.chunk.getUint16();t===18761?this.le=!0:t===19789&&(this.le=!1),this.chunk.le=this.le,this.headerParsed=!0}parseTags(t,e,i=new Map){let{pick:r,skip:s}=this.options[e];r=new Set(r);let o=r.size>0,a=s.size===0,l=this.chunk.getUint16(t);t+=2;for(let u=0;u<l;u++){let h=this.chunk.getUint16(t);if(o){if(r.has(h)&&(i.set(h,this.parseTag(t,h,e)),r.delete(h),r.size===0))break}else!a&&s.has(h)||i.set(h,this.parseTag(t,h,e));t+=12}return i}parseTag(t,e,i){let{chunk:r}=this,s=r.getUint16(t+2),o=r.getUint32(t+4),a=za[s];if(a*o<=4?t+=8:t=r.getUint32(t+8),(s<1||s>13)&&Pt(`Invalid TIFF value type. block: ${i.toUpperCase()}, tag: ${e.toString(16)}, type: ${s}, offset ${t}`),t>r.byteLength&&Pt(`Invalid TIFF value offset. block: ${i.toUpperCase()}, tag: ${e.toString(16)}, type: ${s}, offset ${t} is outside of chunk size ${r.byteLength}`),s===1)return r.getUint8Array(t,o);if(s===2)return Je(r.getString(t,o));if(s===7)return r.getUint8Array(t,o);if(o===1)return this.parseTagValue(s,t);{let l=new(function(h){switch(h){case 1:return Uint8Array;case 3:return Uint16Array;case 4:return Uint32Array;case 5:return Array;case 6:return Int8Array;case 8:return Int16Array;case 9:return Int32Array;case 10:return Array;case 11:return Float32Array;case 12:return Float64Array;default:return Array}}(s))(o),u=a;for(let h=0;h<o;h++)l[h]=this.parseTagValue(s,t),t+=u;return l}}parseTagValue(t,e){let{chunk:i}=this;switch(t){case 1:return i.getUint8(e);case 3:return i.getUint16(e);case 4:return i.getUint32(e);case 5:return i.getUint32(e)/i.getUint32(e+4);case 6:return i.getInt8(e);case 8:return i.getInt16(e);case 9:return i.getInt32(e);case 10:return i.getInt32(e)/i.getInt32(e+4);case 11:return i.getFloat(e);case 12:return i.getDouble(e);case 13:return i.getUint32(e);default:Pt(`Invalid tiff type ${t}`)}}}class Bi extends Va{static canHandle(t,e){return t.getUint8(e+1)===225&&t.getUint32(e+4)===1165519206&&t.getUint16(e+8)===0}async parse(){this.parseHeader();let{options:t}=this;return t.ifd0.enabled&&await this.parseIfd0Block(),t.exif.enabled&&await this.safeParse("parseExifBlock"),t.gps.enabled&&await this.safeParse("parseGpsBlock"),t.interop.enabled&&await this.safeParse("parseInteropBlock"),t.ifd1.enabled&&await this.safeParse("parseThumbnailBlock"),this.createOutput()}safeParse(t){let e=this[t]();return e.catch!==void 0&&(e=e.catch(this.handleError)),e}findIfd0Offset(){this.ifd0Offset===void 0&&(this.ifd0Offset=this.chunk.getUint32(4))}findIfd1Offset(){if(this.ifd1Offset===void 0){this.findIfd0Offset();let t=this.chunk.getUint16(this.ifd0Offset),e=this.ifd0Offset+2+12*t;this.ifd1Offset=this.chunk.getUint32(e)}}parseBlock(t,e){let i=new Map;return this[e]=i,this.parseTags(t,e,i),i}async parseIfd0Block(){if(this.ifd0)return;let{file:t}=this;this.findIfd0Offset(),this.ifd0Offset<8&&Pt("Malformed EXIF data"),!t.chunked&&this.ifd0Offset>t.byteLength&&Pt(`IFD0 offset points to outside of file.
this.ifd0Offset: ${this.ifd0Offset}, file.byteLength: ${t.byteLength}`),t.tiff&&await t.ensureChunk(this.ifd0Offset,Hi(this.options));let e=this.parseBlock(this.ifd0Offset,"ifd0");return e.size!==0?(this.exifOffset=e.get(34665),this.interopOffset=e.get(40965),this.gpsOffset=e.get(34853),this.xmp=e.get(700),this.iptc=e.get(33723),this.icc=e.get(34675),this.options.sanitize&&(e.delete(34665),e.delete(40965),e.delete(34853),e.delete(700),e.delete(33723),e.delete(34675)),e):void 0}async parseExifBlock(){if(this.exif||(this.ifd0||await this.parseIfd0Block(),this.exifOffset===void 0))return;this.file.tiff&&await this.file.ensureChunk(this.exifOffset,Hi(this.options));let t=this.parseBlock(this.exifOffset,"exif");return this.interopOffset||(this.interopOffset=t.get(40965)),this.makerNote=t.get(37500),this.userComment=t.get(37510),this.options.sanitize&&(t.delete(40965),t.delete(37500),t.delete(37510)),this.unpack(t,41728),this.unpack(t,41729),t}unpack(t,e){let i=t.get(e);i&&i.length===1&&t.set(e,i[0])}async parseGpsBlock(){if(this.gps||(this.ifd0||await this.parseIfd0Block(),this.gpsOffset===void 0))return;let t=this.parseBlock(this.gpsOffset,"gps");return t&&t.has(2)&&t.has(4)&&(t.set("latitude",zr(...t.get(2),t.get(1))),t.set("longitude",zr(...t.get(4),t.get(3)))),t}async parseInteropBlock(){if(!this.interop&&(this.ifd0||await this.parseIfd0Block(),this.interopOffset!==void 0||this.exif||await this.parseExifBlock(),this.interopOffset!==void 0))return this.parseBlock(this.interopOffset,"interop")}async parseThumbnailBlock(t=!1){if(!this.ifd1&&!this.ifd1Parsed&&(!this.options.mergeOutput||t))return this.findIfd1Offset(),this.ifd1Offset>0&&(this.parseBlock(this.ifd1Offset,"ifd1"),this.ifd1Parsed=!0),this.ifd1}async extractThumbnail(){if(this.headerParsed||this.parseHeader(),this.ifd1Parsed||await this.parseThumbnailBlock(!0),this.ifd1===void 0)return;let t=this.ifd1.get(513),e=this.ifd1.get(514);return this.chunk.getUint8Array(t,e)}get image(){return this.ifd0}get thumbnail(){return this.ifd1}createOutput(){let t,e,i,r={};for(e of bt)if(t=this[e],!gs(t))if(i=this.canTranslate?this.translateBlock(t,e):Object.fromEntries(t),this.options.mergeOutput){if(e==="ifd1")continue;Object.assign(r,i)}else r[e]=i;return this.makerNote&&(r.makerNote=this.makerNote),this.userComment&&(r.userComment=this.userComment),r}assignToOutput(t,e){if(this.globalOptions.mergeOutput)Object.assign(t,e);else for(let[i,r]of Object.entries(e))this.assignObjectToOutput(t,i,r)}}function zr(n,t,e,i){var r=n+t/60+e/3600;return i!=="S"&&i!=="W"||(r*=-1),r}tt(Bi,"type","tiff"),tt(Bi,"headerLength",10),Mt.set("tiff",Bi);var Ga=Object.freeze({__proto__:null,default:Ua,Exifr:Re,fileParsers:he,segmentParsers:Mt,fileReaders:fe,tagKeys:kt,tagValues:$t,tagRevivers:Ne,createDictionary:wt,extendDictionary:on,fetchUrlAsArrayBuffer:sn,readBlobAsArrayBuffer:an,chunkedProps:Ce,otherSegments:yn,segments:ln,tiffBlocks:bt,segmentsAndBlocks:Fe,tiffExtractables:ke,inheritables:xn,allFormatters:Ae,Options:cn,parse:gi});const ir={ifd0:!1,ifd1:!1,exif:!1,gps:!1,interop:!1,sanitize:!1,reviveValues:!0,translateKeys:!1,translateValues:!1,mergeOutput:!1},rr=Object.assign({},ir,{firstChunkSize:4e4,gps:[1,2,3,4]});async function bs(n){let t=new Re(rr);await t.read(n);let e=await t.parse();if(e&&e.gps){let{latitude:i,longitude:r}=e.gps;return{latitude:i,longitude:r}}}const sr=Object.assign({},ir,{tiff:!1,ifd1:!0,mergeOutput:!1});async function _s(n){let t=new Re(sr);await t.read(n);let e=await t.extractThumbnail();return e&&mi?pi.from(e):e}async function ws(n){let t=await this.thumbnail(n);if(t!==void 0){let e=new Blob([t]);return URL.createObjectURL(e)}}const ar=Object.assign({},ir,{firstChunkSize:4e4,ifd0:[274]});async function or(n){let t=new Re(ar);await t.read(n);let e=await t.parse();if(e&&e.ifd0)return e.ifd0[274]}const lr=Object.freeze({1:{dimensionSwapped:!1,scaleX:1,scaleY:1,deg:0,rad:0},2:{dimensionSwapped:!1,scaleX:-1,scaleY:1,deg:0,rad:0},3:{dimensionSwapped:!1,scaleX:1,scaleY:1,deg:180,rad:180*Math.PI/180},4:{dimensionSwapped:!1,scaleX:-1,scaleY:1,deg:180,rad:180*Math.PI/180},5:{dimensionSwapped:!0,scaleX:1,scaleY:-1,deg:90,rad:90*Math.PI/180},6:{dimensionSwapped:!0,scaleX:1,scaleY:1,deg:90,rad:90*Math.PI/180},7:{dimensionSwapped:!0,scaleX:1,scaleY:-1,deg:270,rad:270*Math.PI/180},8:{dimensionSwapped:!0,scaleX:1,scaleY:1,deg:270,rad:270*Math.PI/180}});let Ze=!0,tn=!0;if(typeof navigator=="object"){let n=navigator.userAgent;if(n.includes("iPad")||n.includes("iPhone")){let t=n.match(/OS (\d+)_(\d+)/);if(t){let[,e,i]=t;Ze=Number(e)+.1*Number(i)<13.4,tn=!1}}else if(n.includes("OS X 10")){let[,t]=n.match(/OS X 10[_.](\d+)/);Ze=tn=Number(t)<15}if(n.includes("Chrome/")){let[,t]=n.match(/Chrome\/(\d+)/);Ze=tn=Number(t)<81}else if(n.includes("Firefox/")){let[,t]=n.match(/Firefox\/(\d+)/);Ze=tn=Number(t)<77}}async function Ms(n){let t=await or(n);return Object.assign({canvas:Ze,css:tn},lr[t])}class Xa extends zt{constructor(...t){super(...t),tt(this,"ranges",new Ya),this.byteLength!==0&&this.ranges.add(0,this.byteLength)}_tryExtend(t,e,i){if(t===0&&this.byteLength===0&&i){let r=new DataView(i.buffer||i,i.byteOffset,i.byteLength);this._swapDataView(r)}else{let r=t+e;if(r>this.byteLength){let{dataView:s}=this._extend(r);this._swapDataView(s)}}}_extend(t){let e;e=mi?pi.allocUnsafe(t):new Uint8Array(t);let i=new DataView(e.buffer,e.byteOffset,e.byteLength);return e.set(new Uint8Array(this.buffer,this.byteOffset,this.byteLength),0),{uintView:e,dataView:i}}subarray(t,e,i=!1){return e=e||this._lengthToEnd(t),i&&this._tryExtend(t,e),this.ranges.add(t,e),super.subarray(t,e)}set(t,e,i=!1){i&&this._tryExtend(e,t.byteLength,t);let r=super.set(t,e);return this.ranges.add(e,r.byteLength),r}async ensureChunk(t,e){this.chunked&&(this.ranges.available(t,e)||await this.readChunk(t,e))}available(t,e){return this.ranges.available(t,e)}}class Ya{constructor(){tt(this,"list",[])}get length(){return this.list.length}add(t,e,i=0){let r=t+e,s=this.list.filter(o=>Vr(t,o.offset,r)||Vr(t,o.end,r));if(s.length>0){t=Math.min(t,...s.map(a=>a.offset)),r=Math.max(r,...s.map(a=>a.end)),e=r-t;let o=s.shift();o.offset=t,o.length=e,o.end=r,this.list=this.list.filter(a=>!s.includes(a))}else this.list.push({offset:t,length:e,end:r})}available(t,e){let i=t+e;return this.list.some(r=>r.offset<=t&&i<=r.end)}}function Vr(n,t,e){return n<=t&&t<=e}class xi extends Xa{constructor(t,e){super(0),tt(this,"chunksRead",0),this.input=t,this.options=e}async readWhole(){this.chunked=!1,await this.readChunk(this.nextChunkOffset)}async readChunked(){this.chunked=!0,await this.readChunk(0,this.options.firstChunkSize)}async readNextChunk(t=this.nextChunkOffset){if(this.fullyRead)return this.chunksRead++,!1;let e=this.options.chunkSize,i=await this.readChunk(t,e);return!!i&&i.byteLength===e}async readChunk(t,e){if(this.chunksRead++,(e=this.safeWrapAddress(t,e))!==0)return this._readChunk(t,e)}safeWrapAddress(t,e){return this.size!==void 0&&t+e>this.size?Math.max(0,this.size-t):e}get nextChunkOffset(){if(this.ranges.list.length!==0)return this.ranges.list[0].length}get canReadNextChunk(){return this.chunksRead<this.options.chunkLimit}get fullyRead(){return this.size!==void 0&&this.nextChunkOffset===this.size}read(){return this.options.chunked?this.readChunked():this.readWhole()}close(){}}fe.set("blob",class extends xi{async readWhole(){this.chunked=!1;let n=await an(this.input);this._swapArrayBuffer(n)}readChunked(){return this.chunked=!0,this.size=this.input.size,super.readChunked()}async _readChunk(n,t){let e=t?n+t:void 0,i=this.input.slice(n,e),r=await an(i);return this.set(r,n,!0)}});var Wa=Object.freeze({__proto__:null,default:Ga,Exifr:Re,fileParsers:he,segmentParsers:Mt,fileReaders:fe,tagKeys:kt,tagValues:$t,tagRevivers:Ne,createDictionary:wt,extendDictionary:on,fetchUrlAsArrayBuffer:sn,readBlobAsArrayBuffer:an,chunkedProps:Ce,otherSegments:yn,segments:ln,tiffBlocks:bt,segmentsAndBlocks:Fe,tiffExtractables:ke,inheritables:xn,allFormatters:Ae,Options:cn,parse:gi,gpsOnlyOptions:rr,gps:bs,thumbnailOnlyOptions:sr,thumbnail:_s,thumbnailUrl:ws,orientationOnlyOptions:ar,orientation:or,rotations:lr,get rotateCanvas(){return Ze},get rotateCss(){return tn},rotation:Ms});fe.set("url",class extends xi{async readWhole(){this.chunked=!1;let n=await sn(this.input);n instanceof ArrayBuffer?this._swapArrayBuffer(n):n instanceof Uint8Array&&this._swapBuffer(n)}async _readChunk(n,t){let e=t?n+t-1:void 0,i=this.options.httpHeaders||{};(n||e)&&(i.range=`bytes=${[n,e].join("-")}`);let r=await er(this.input,{headers:i}),s=await r.arrayBuffer(),o=s.byteLength;if(r.status!==416)return o!==t&&(this.size=n+o),this.set(s,n,!0)}});zt.prototype.getUint64=function(n){let t=this.getUint32(n),e=this.getUint32(n+4);return t<1048575?t<<32|e:typeof Hn!==void 0?(console.warn("Using BigInt because of type 64uint but JS can only handle 53b numbers."),Hn(t)<<Hn(32)|Hn(e)):void Pt("Trying to read 64b value but JS can only handle 53b numbers.")};class Ha extends yi{parseBoxes(t=0){let e=[];for(;t<this.file.byteLength-4;){let i=this.parseBoxHead(t);if(e.push(i),i.length===0)break;t+=i.length}return e}parseSubBoxes(t){t.boxes=this.parseBoxes(t.start)}findBox(t,e){return t.boxes===void 0&&this.parseSubBoxes(t),t.boxes.find(i=>i.kind===e)}parseBoxHead(t){let e=this.file.getUint32(t),i=this.file.getString(t+4,4),r=t+8;return e===1&&(e=this.file.getUint64(t+8),r+=8),{offset:t,length:e,kind:i,start:r}}parseBoxFullHead(t){if(t.version!==void 0)return;let e=this.file.getUint32(t.start);t.version=e>>24,t.start+=4}}class Ss extends Ha{static canHandle(t,e){if(e!==0)return!1;let i=t.getUint16(2);if(i>50)return!1;let r=16,s=[];for(;r<i;)s.push(t.getString(r,4)),r+=4;return s.includes(this.type)}async parse(){let t=this.file.getUint32(0),e=this.parseBoxHead(t);for(;e.kind!=="meta";)t+=e.length,await this.file.ensureChunk(t,16),e=this.parseBoxHead(t);await this.file.ensureChunk(e.offset,e.length),this.parseBoxFullHead(e),this.parseSubBoxes(e),this.options.icc.enabled&&await this.findIcc(e),this.options.tiff.enabled&&await this.findExif(e)}async registerSegment(t,e,i){await this.file.ensureChunk(e,i);let r=this.file.subarray(e,i);this.createParser(t,r)}async findIcc(t){let e=this.findBox(t,"iprp");if(e===void 0)return;let i=this.findBox(e,"ipco");if(i===void 0)return;let r=this.findBox(i,"colr");r!==void 0&&await this.registerSegment("icc",r.offset+12,r.length)}async findExif(t){let e=this.findBox(t,"iinf");if(e===void 0)return;let i=this.findBox(t,"iloc");if(i===void 0)return;let r=this.findExifLocIdInIinf(e),s=this.findExtentInIloc(i,r);if(s===void 0)return;let[o,a]=s;await this.file.ensureChunk(o,a);let l=4+this.file.getUint32(o);o+=l,a-=l,await this.registerSegment("tiff",o,a)}findExifLocIdInIinf(t){this.parseBoxFullHead(t);let e,i,r,s,o=t.start,a=this.file.getUint16(o);for(o+=2;a--;){if(e=this.parseBoxHead(o),this.parseBoxFullHead(e),i=e.start,e.version>=2&&(r=e.version===3?4:2,s=this.file.getString(i+r+2,4),s==="Exif"))return this.file.getUintBytes(i,r);o+=e.length}}get8bits(t){let e=this.file.getUint8(t);return[e>>4,15&e]}findExtentInIloc(t,e){this.parseBoxFullHead(t);let i=t.start,[r,s]=this.get8bits(i++),[o,a]=this.get8bits(i++),l=t.version===2?4:2,u=t.version===1||t.version===2?2:0,h=a+r+s,c=t.version===2?4:2,m=this.file.getUintBytes(i,c);for(i+=c;m--;){let d=this.file.getUintBytes(i,l);i+=l+u+2+o;let p=this.file.getUint16(i);if(i+=2,d===e)return p>1&&console.warn(`ILOC box has more than one extent but we're only processing one
Please create an issue at https://github.com/MikeKovarik/exifr with this file`),[this.file.getUintBytes(i+a,r),this.file.getUintBytes(i+a+r,s)];i+=p*h}}}class Ps extends Ss{}tt(Ps,"type","heic");class Gr extends Ss{}tt(Gr,"type","avif"),he.set("heic",Ps),he.set("avif",Gr),wt(kt,["ifd0","ifd1"],[[256,"ImageWidth"],[257,"ImageHeight"],[258,"BitsPerSample"],[259,"Compression"],[262,"PhotometricInterpretation"],[270,"ImageDescription"],[271,"Make"],[272,"Model"],[273,"StripOffsets"],[274,"Orientation"],[277,"SamplesPerPixel"],[278,"RowsPerStrip"],[279,"StripByteCounts"],[282,"XResolution"],[283,"YResolution"],[284,"PlanarConfiguration"],[296,"ResolutionUnit"],[301,"TransferFunction"],[305,"Software"],[306,"ModifyDate"],[315,"Artist"],[316,"HostComputer"],[317,"Predictor"],[318,"WhitePoint"],[319,"PrimaryChromaticities"],[513,"ThumbnailOffset"],[514,"ThumbnailLength"],[529,"YCbCrCoefficients"],[530,"YCbCrSubSampling"],[531,"YCbCrPositioning"],[532,"ReferenceBlackWhite"],[700,"ApplicationNotes"],[33432,"Copyright"],[33723,"IPTC"],[34665,"ExifIFD"],[34675,"ICC"],[34853,"GpsIFD"],[330,"SubIFD"],[40965,"InteropIFD"],[40091,"XPTitle"],[40092,"XPComment"],[40093,"XPAuthor"],[40094,"XPKeywords"],[40095,"XPSubject"]]),wt(kt,"exif",[[33434,"ExposureTime"],[33437,"FNumber"],[34850,"ExposureProgram"],[34852,"SpectralSensitivity"],[34855,"ISO"],[34858,"TimeZoneOffset"],[34859,"SelfTimerMode"],[34864,"SensitivityType"],[34865,"StandardOutputSensitivity"],[34866,"RecommendedExposureIndex"],[34867,"ISOSpeed"],[34868,"ISOSpeedLatitudeyyy"],[34869,"ISOSpeedLatitudezzz"],[36864,"ExifVersion"],[36867,"DateTimeOriginal"],[36868,"CreateDate"],[36873,"GooglePlusUploadCode"],[36880,"OffsetTime"],[36881,"OffsetTimeOriginal"],[36882,"OffsetTimeDigitized"],[37121,"ComponentsConfiguration"],[37122,"CompressedBitsPerPixel"],[37377,"ShutterSpeedValue"],[37378,"ApertureValue"],[37379,"BrightnessValue"],[37380,"ExposureCompensation"],[37381,"MaxApertureValue"],[37382,"SubjectDistance"],[37383,"MeteringMode"],[37384,"LightSource"],[37385,"Flash"],[37386,"FocalLength"],[37393,"ImageNumber"],[37394,"SecurityClassification"],[37395,"ImageHistory"],[37396,"SubjectArea"],[37500,"MakerNote"],[37510,"UserComment"],[37520,"SubSecTime"],[37521,"SubSecTimeOriginal"],[37522,"SubSecTimeDigitized"],[37888,"AmbientTemperature"],[37889,"Humidity"],[37890,"Pressure"],[37891,"WaterDepth"],[37892,"Acceleration"],[37893,"CameraElevationAngle"],[40960,"FlashpixVersion"],[40961,"ColorSpace"],[40962,"ExifImageWidth"],[40963,"ExifImageHeight"],[40964,"RelatedSoundFile"],[41483,"FlashEnergy"],[41486,"FocalPlaneXResolution"],[41487,"FocalPlaneYResolution"],[41488,"FocalPlaneResolutionUnit"],[41492,"SubjectLocation"],[41493,"ExposureIndex"],[41495,"SensingMethod"],[41728,"FileSource"],[41729,"SceneType"],[41730,"CFAPattern"],[41985,"CustomRendered"],[41986,"ExposureMode"],[41987,"WhiteBalance"],[41988,"DigitalZoomRatio"],[41989,"FocalLengthIn35mmFormat"],[41990,"SceneCaptureType"],[41991,"GainControl"],[41992,"Contrast"],[41993,"Saturation"],[41994,"Sharpness"],[41996,"SubjectDistanceRange"],[42016,"ImageUniqueID"],[42032,"OwnerName"],[42033,"SerialNumber"],[42034,"LensInfo"],[42035,"LensMake"],[42036,"LensModel"],[42037,"LensSerialNumber"],[42080,"CompositeImage"],[42081,"CompositeImageCount"],[42082,"CompositeImageExposureTimes"],[42240,"Gamma"],[59932,"Padding"],[59933,"OffsetSchema"],[65e3,"OwnerName"],[65001,"SerialNumber"],[65002,"Lens"],[65100,"RawFile"],[65101,"Converter"],[65102,"WhiteBalance"],[65105,"Exposure"],[65106,"Shadows"],[65107,"Brightness"],[65108,"Contrast"],[65109,"Saturation"],[65110,"Sharpness"],[65111,"Smoothness"],[65112,"MoireFilter"],[40965,"InteropIFD"]]),wt(kt,"gps",[[0,"GPSVersionID"],[1,"GPSLatitudeRef"],[2,"GPSLatitude"],[3,"GPSLongitudeRef"],[4,"GPSLongitude"],[5,"GPSAltitudeRef"],[6,"GPSAltitude"],[7,"GPSTimeStamp"],[8,"GPSSatellites"],[9,"GPSStatus"],[10,"GPSMeasureMode"],[11,"GPSDOP"],[12,"GPSSpeedRef"],[13,"GPSSpeed"],[14,"GPSTrackRef"],[15,"GPSTrack"],[16,"GPSImgDirectionRef"],[17,"GPSImgDirection"],[18,"GPSMapDatum"],[19,"GPSDestLatitudeRef"],[20,"GPSDestLatitude"],[21,"GPSDestLongitudeRef"],[22,"GPSDestLongitude"],[23,"GPSDestBearingRef"],[24,"GPSDestBearing"],[25,"GPSDestDistanceRef"],[26,"GPSDestDistance"],[27,"GPSProcessingMethod"],[28,"GPSAreaInformation"],[29,"GPSDateStamp"],[30,"GPSDifferential"],[31,"GPSHPositioningError"]]),wt($t,["ifd0","ifd1"],[[274,{1:"Horizontal (normal)",2:"Mirror horizontal",3:"Rotate 180",4:"Mirror vertical",5:"Mirror horizontal and rotate 270 CW",6:"Rotate 90 CW",7:"Mirror horizontal and rotate 90 CW",8:"Rotate 270 CW"}],[296,{1:"None",2:"inches",3:"cm"}]]);let In=wt($t,"exif",[[34850,{0:"Not defined",1:"Manual",2:"Normal program",3:"Aperture priority",4:"Shutter priority",5:"Creative program",6:"Action program",7:"Portrait mode",8:"Landscape mode"}],[37121,{0:"-",1:"Y",2:"Cb",3:"Cr",4:"R",5:"G",6:"B"}],[37383,{0:"Unknown",1:"Average",2:"CenterWeightedAverage",3:"Spot",4:"MultiSpot",5:"Pattern",6:"Partial",255:"Other"}],[37384,{0:"Unknown",1:"Daylight",2:"Fluorescent",3:"Tungsten (incandescent light)",4:"Flash",9:"Fine weather",10:"Cloudy weather",11:"Shade",12:"Daylight fluorescent (D 5700 - 7100K)",13:"Day white fluorescent (N 4600 - 5400K)",14:"Cool white fluorescent (W 3900 - 4500K)",15:"White fluorescent (WW 3200 - 3700K)",17:"Standard light A",18:"Standard light B",19:"Standard light C",20:"D55",21:"D65",22:"D75",23:"D50",24:"ISO studio tungsten",255:"Other"}],[37385,{0:"Flash did not fire",1:"Flash fired",5:"Strobe return light not detected",7:"Strobe return light detected",9:"Flash fired, compulsory flash mode",13:"Flash fired, compulsory flash mode, return light not detected",15:"Flash fired, compulsory flash mode, return light detected",16:"Flash did not fire, compulsory flash mode",24:"Flash did not fire, auto mode",25:"Flash fired, auto mode",29:"Flash fired, auto mode, return light not detected",31:"Flash fired, auto mode, return light detected",32:"No flash function",65:"Flash fired, red-eye reduction mode",69:"Flash fired, red-eye reduction mode, return light not detected",71:"Flash fired, red-eye reduction mode, return light detected",73:"Flash fired, compulsory flash mode, red-eye reduction mode",77:"Flash fired, compulsory flash mode, red-eye reduction mode, return light not detected",79:"Flash fired, compulsory flash mode, red-eye reduction mode, return light detected",89:"Flash fired, auto mode, red-eye reduction mode",93:"Flash fired, auto mode, return light not detected, red-eye reduction mode",95:"Flash fired, auto mode, return light detected, red-eye reduction mode"}],[41495,{1:"Not defined",2:"One-chip color area sensor",3:"Two-chip color area sensor",4:"Three-chip color area sensor",5:"Color sequential area sensor",7:"Trilinear sensor",8:"Color sequential linear sensor"}],[41728,{1:"Film Scanner",2:"Reflection Print Scanner",3:"Digital Camera"}],[41729,{1:"Directly photographed"}],[41985,{0:"Normal",1:"Custom",2:"HDR (no original saved)",3:"HDR (original saved)",4:"Original (for HDR)",6:"Panorama",7:"Portrait HDR",8:"Portrait"}],[41986,{0:"Auto",1:"Manual",2:"Auto bracket"}],[41987,{0:"Auto",1:"Manual"}],[41990,{0:"Standard",1:"Landscape",2:"Portrait",3:"Night",4:"Other"}],[41991,{0:"None",1:"Low gain up",2:"High gain up",3:"Low gain down",4:"High gain down"}],[41996,{0:"Unknown",1:"Macro",2:"Close",3:"Distant"}],[42080,{0:"Unknown",1:"Not a Composite Image",2:"General Composite Image",3:"Composite Image Captured While Shooting"}]]);const Xr={1:"No absolute unit of measurement",2:"Inch",3:"Centimeter"};In.set(37392,Xr),In.set(41488,Xr);const Oi={0:"Normal",1:"Low",2:"High"};function Yr(n){return typeof n=="object"&&n.length!==void 0?n[0]:n}function Wr(n){let t=Array.from(n).slice(1);return t[1]>15&&(t=t.map(e=>String.fromCharCode(e))),t[2]!=="0"&&t[2]!==0||t.pop(),t.join(".")}function zi(n){if(typeof n=="string"){var[t,e,i,r,s,o]=n.trim().split(/[-: ]/g).map(Number),a=new Date(t,e-1,i);return Number.isNaN(r)||Number.isNaN(s)||Number.isNaN(o)||(a.setHours(r),a.setMinutes(s),a.setSeconds(o)),Number.isNaN(+a)?n:a}}function kn(n){if(typeof n=="string")return n;let t=[];if(n[1]===0&&n[n.length-1]===0)for(let e=0;e<n.length;e+=2)t.push(Hr(n[e+1],n[e]));else for(let e=0;e<n.length;e+=2)t.push(Hr(n[e],n[e+1]));return Je(String.fromCodePoint(...t))}function Hr(n,t){return n<<8|t}In.set(41992,Oi),In.set(41993,Oi),In.set(41994,Oi),wt(Ne,["ifd0","ifd1"],[[50827,function(n){return typeof n!="string"?ys(n):n}],[306,zi],[40091,kn],[40092,kn],[40093,kn],[40094,kn],[40095,kn]]),wt(Ne,"exif",[[40960,Wr],[36864,Wr],[36867,zi],[36868,zi],[40962,Yr],[40963,Yr]]),wt(Ne,"gps",[[0,n=>Array.from(n).join(".")],[7,n=>Array.from(n).join(":")]]);class Vi extends Kt{static canHandle(t,e){return t.getUint8(e+1)===225&&t.getUint32(e+4)===1752462448&&t.getString(e+4,20)==="http://ns.adobe.com/"}static headerLength(t,e){return t.getString(e+4,34)==="http://ns.adobe.com/xmp/extension/"?79:33}static findPosition(t,e){let i=super.findPosition(t,e);return i.multiSegment=i.extended=i.headerLength===79,i.multiSegment?(i.chunkCount=t.getUint8(e+72),i.chunkNumber=t.getUint8(e+76),t.getUint8(e+77)!==0&&i.chunkNumber++):(i.chunkCount=1/0,i.chunkNumber=-1),i}static handleMultiSegments(t){return t.map(e=>e.chunk.getString()).join("")}normalizeInput(t){return typeof t=="string"?t:zt.from(t).getString()}parse(t=this.chunk){if(!this.localOptions.parse)return t;t=function(s){let o={},a={};for(let l of ks)o[l]=[],a[l]=0;return s.replace(Ka,(l,u,h)=>{if(u==="<"){let c=++a[h];return o[h].push(c),`${l}#${c}`}return`${l}#${o[h].pop()}`})}(t);let e=pn.findAll(t,"rdf","Description");e.length===0&&e.push(new pn("rdf","Description",void 0,t));let i,r={};for(let s of e)for(let o of s.properties)i=Qa(o.ns,r),vs(o,i);return function(s){let o;for(let a in s)o=s[a]=ii(s[a]),o===void 0&&delete s[a];return ii(s)}(r)}assignToOutput(t,e){if(this.localOptions.parse)for(let[i,r]of Object.entries(e))switch(i){case"tiff":this.assignObjectToOutput(t,"ifd0",r);break;case"exif":this.assignObjectToOutput(t,"exif",r);break;case"xmlns":break;default:this.assignObjectToOutput(t,i,r)}else t.xmp=e}}tt(Vi,"type","xmp"),tt(Vi,"multiSegment",!0),Mt.set("xmp",Vi);class si{static findAll(t){return Cs(t,/([a-zA-Z0-9-]+):([a-zA-Z0-9-]+)=("[^"]*"|'[^']*')/gm).map(si.unpackMatch)}static unpackMatch(t){let e=t[1],i=t[2],r=t[3].slice(1,-1);return r=Fs(r),new si(e,i,r)}constructor(t,e,i){this.ns=t,this.name=e,this.value=i}serialize(){return this.value}}class pn{static findAll(t,e,i){if(e!==void 0||i!==void 0){e=e||"[\\w\\d-]+",i=i||"[\\w\\d-]+";var r=new RegExp(`<(${e}):(${i})(#\\d+)?((\\s+?[\\w\\d-:]+=("[^"]*"|'[^']*'))*\\s*)(\\/>|>([\\s\\S]*?)<\\/\\1:\\2\\3>)`,"gm")}else r=/<([\w\d-]+):([\w\d-]+)(#\d+)?((\s+?[\w\d-:]+=("[^"]*"|'[^']*'))*\s*)(\/>|>([\s\S]*?)<\/\1:\2\3>)/gm;return Cs(t,r).map(pn.unpackMatch)}static unpackMatch(t){let e=t[1],i=t[2],r=t[4],s=t[8];return new pn(e,i,r,s)}constructor(t,e,i,r){this.ns=t,this.name=e,this.attrString=i,this.innerXml=r,this.attrs=si.findAll(i),this.children=pn.findAll(r),this.value=this.children.length===0?Fs(r):void 0,this.properties=[...this.attrs,...this.children]}get isPrimitive(){return this.value!==void 0&&this.attrs.length===0&&this.children.length===0}get isListContainer(){return this.children.length===1&&this.children[0].isList}get isList(){let{ns:t,name:e}=this;return t==="rdf"&&(e==="Seq"||e==="Bag"||e==="Alt")}get isListItem(){return this.ns==="rdf"&&this.name==="li"}serialize(){if(this.properties.length===0&&this.value===void 0)return;if(this.isPrimitive)return this.value;if(this.isListContainer)return this.children[0].serialize();if(this.isList)return qa(this.children.map(ja));if(this.isListItem&&this.children.length===1&&this.attrs.length===0)return this.children[0].serialize();let t={};for(let e of this.properties)vs(e,t);return this.value!==void 0&&(t.value=this.value),ii(t)}}function vs(n,t){let e=n.serialize();e!==void 0&&(t[n.name]=e)}var ja=n=>n.serialize(),qa=n=>n.length===1?n[0]:n,Qa=(n,t)=>t[n]?t[n]:t[n]={};function Cs(n,t){let e,i=[];if(!n)return i;for(;(e=t.exec(n))!==null;)i.push(e);return i}function Fs(n){if(function(i){return i==null||i==="null"||i==="undefined"||i===""||i.trim()===""}(n))return;let t=Number(n);if(!Number.isNaN(t))return t;let e=n.toLowerCase();return e==="true"||e!=="false"&&n.trim()}const ks=["rdf:li","rdf:Seq","rdf:Bag","rdf:Alt","rdf:Description"],Ka=new RegExp(`(<|\\/)(${ks.join("|")})`,"g");var As=Object.freeze({__proto__:null,default:Wa,Exifr:Re,fileParsers:he,segmentParsers:Mt,fileReaders:fe,tagKeys:kt,tagValues:$t,tagRevivers:Ne,createDictionary:wt,extendDictionary:on,fetchUrlAsArrayBuffer:sn,readBlobAsArrayBuffer:an,chunkedProps:Ce,otherSegments:yn,segments:ln,tiffBlocks:bt,segmentsAndBlocks:Fe,tiffExtractables:ke,inheritables:xn,allFormatters:Ae,Options:cn,parse:gi,gpsOnlyOptions:rr,gps:bs,thumbnailOnlyOptions:sr,thumbnail:_s,thumbnailUrl:ws,orientationOnlyOptions:ar,orientation:or,rotations:lr,get rotateCanvas(){return Ze},get rotateCss(){return tn},rotation:Ms});let jr=ni("fs",n=>n.promises);fe.set("fs",class extends xi{async readWhole(){this.chunked=!1,this.fs=await jr;let n=await this.fs.readFile(this.input);this._swapBuffer(n)}async readChunked(){this.chunked=!0,this.fs=await jr,await this.open(),await this.readChunk(0,this.options.firstChunkSize)}async open(){this.fh===void 0&&(this.fh=await this.fs.open(this.input,"r"),this.size=(await this.fh.stat(this.input)).size)}async _readChunk(n,t){this.fh===void 0&&await this.open(),n+t>this.size&&(t=this.size-n);var e=this.subarray(n,t,!0);return await this.fh.read(e.dataView,0,t,n),e}async close(){if(this.fh){let n=this.fh;this.fh=void 0,await n.close()}}});fe.set("base64",class extends xi{constructor(...n){super(...n),this.input=this.input.replace(/^data:([^;]+);base64,/gim,""),this.size=this.input.length/4*3,this.input.endsWith("==")?this.size-=2:this.input.endsWith("=")&&(this.size-=1)}async _readChunk(n,t){let e,i,r=this.input;n===void 0?(n=0,e=0,i=0):(e=4*Math.floor(n/3),i=n-e/4*3),t===void 0&&(t=this.size);let s=n+t,o=e+4*Math.ceil(s/3);r=r.slice(e,o);let a=Math.min(t,this.size-n);if(mi){let l=pi.from(r,"base64").slice(i,i+a);return this.set(l,n,!0)}{let l=this.subarray(n,a,!0),u=atob(r),h=l.toUint8();for(let c=0;c<a;c++)h[c]=u.charCodeAt(i+c);return l}}});class qr extends yi{static canHandle(t,e){return e===18761||e===19789}extendOptions(t){let{ifd0:e,xmp:i,iptc:r,icc:s}=t;i.enabled&&e.deps.add(700),r.enabled&&e.deps.add(33723),s.enabled&&e.deps.add(34675),e.finalizeFilters()}async parse(){let{tiff:t,xmp:e,iptc:i,icc:r}=this.options;if(t.enabled||e.enabled||i.enabled||r.enabled){let s=Math.max(Hi(this.options),this.options.chunkSize);await this.file.ensureChunk(0,s),this.createParser("tiff",this.file),this.parsers.tiff.parseHeader(),await this.parsers.tiff.parseIfd0Block(),this.adaptTiffPropAsSegment("xmp"),this.adaptTiffPropAsSegment("iptc"),this.adaptTiffPropAsSegment("icc")}}adaptTiffPropAsSegment(t){if(this.parsers.tiff[t]){let e=this.parsers.tiff[t];this.injectSegment(t,e)}}}tt(qr,"type","tiff"),he.set("tiff",qr);let $a=ni("zlib");const Ja=["ihdr","iccp","text","itxt","exif"];class Qr extends yi{constructor(...t){super(...t),tt(this,"catchError",e=>this.errors.push(e)),tt(this,"metaChunks",[]),tt(this,"unknownChunks",[])}static canHandle(t,e){return e===35152&&t.getUint32(0)===2303741511&&t.getUint32(4)===218765834}async parse(){let{file:t}=this;await this.findPngChunksInRange(8,t.byteLength),await this.readSegments(this.metaChunks),this.findIhdr(),this.parseTextChunks(),await this.findExif().catch(this.catchError),await this.findXmp().catch(this.catchError),await this.findIcc().catch(this.catchError)}async findPngChunksInRange(t,e){let{file:i}=this;for(;t<e;){let r=i.getUint32(t),s=i.getUint32(t+4),o=i.getString(t+4,4).toLowerCase(),a=r+4+4+4,l={type:o,offset:t,length:a,start:t+4+4,size:r,marker:s};Ja.includes(o)?this.metaChunks.push(l):this.unknownChunks.push(l),t+=a}}parseTextChunks(){let t=this.metaChunks.filter(e=>e.type==="text");for(let e of t){let[i,r]=this.file.getString(e.start,e.size).split("\0");this.injectKeyValToIhdr(i,r)}}injectKeyValToIhdr(t,e){let i=this.parsers.ihdr;i&&i.raw.set(t,e)}findIhdr(){let t=this.metaChunks.find(e=>e.type==="ihdr");t&&this.options.ihdr.enabled!==!1&&this.createParser("ihdr",t.chunk)}async findExif(){let t=this.metaChunks.find(e=>e.type==="exif");t&&this.injectSegment("tiff",t.chunk)}async findXmp(){let t=this.metaChunks.filter(e=>e.type==="itxt");for(let e of t)e.chunk.getString(0,17)==="XML:com.adobe.xmp"&&this.injectSegment("xmp",e.chunk)}async findIcc(){let t=this.metaChunks.find(a=>a.type==="iccp");if(!t)return;let{chunk:e}=t,i=e.getUint8Array(0,81),r=0;for(;r<80&&i[r]!==0;)r++;let s=r+2,o=e.getString(0,r);if(this.injectKeyValToIhdr("ProfileName",o),ei){let a=await $a,l=e.getUint8Array(s);l=a.inflateSync(l),this.injectSegment("icc",l)}}}tt(Qr,"type","png"),he.set("png",Qr),wt(kt,"interop",[[1,"InteropIndex"],[2,"InteropVersion"],[4096,"RelatedImageFileFormat"],[4097,"RelatedImageWidth"],[4098,"RelatedImageHeight"]]),on(kt,"ifd0",[[11,"ProcessingSoftware"],[254,"SubfileType"],[255,"OldSubfileType"],[263,"Thresholding"],[264,"CellWidth"],[265,"CellLength"],[266,"FillOrder"],[269,"DocumentName"],[280,"MinSampleValue"],[281,"MaxSampleValue"],[285,"PageName"],[286,"XPosition"],[287,"YPosition"],[290,"GrayResponseUnit"],[297,"PageNumber"],[321,"HalftoneHints"],[322,"TileWidth"],[323,"TileLength"],[332,"InkSet"],[337,"TargetPrinter"],[18246,"Rating"],[18249,"RatingPercent"],[33550,"PixelScale"],[34264,"ModelTransform"],[34377,"PhotoshopSettings"],[50706,"DNGVersion"],[50707,"DNGBackwardVersion"],[50708,"UniqueCameraModel"],[50709,"LocalizedCameraModel"],[50736,"DNGLensInfo"],[50739,"ShadowScale"],[50740,"DNGPrivateData"],[33920,"IntergraphMatrix"],[33922,"ModelTiePoint"],[34118,"SEMInfo"],[34735,"GeoTiffDirectory"],[34736,"GeoTiffDoubleParams"],[34737,"GeoTiffAsciiParams"],[50341,"PrintIM"],[50721,"ColorMatrix1"],[50722,"ColorMatrix2"],[50723,"CameraCalibration1"],[50724,"CameraCalibration2"],[50725,"ReductionMatrix1"],[50726,"ReductionMatrix2"],[50727,"AnalogBalance"],[50728,"AsShotNeutral"],[50729,"AsShotWhiteXY"],[50730,"BaselineExposure"],[50731,"BaselineNoise"],[50732,"BaselineSharpness"],[50734,"LinearResponseLimit"],[50735,"CameraSerialNumber"],[50741,"MakerNoteSafety"],[50778,"CalibrationIlluminant1"],[50779,"CalibrationIlluminant2"],[50781,"RawDataUniqueID"],[50827,"OriginalRawFileName"],[50828,"OriginalRawFileData"],[50831,"AsShotICCProfile"],[50832,"AsShotPreProfileMatrix"],[50833,"CurrentICCProfile"],[50834,"CurrentPreProfileMatrix"],[50879,"ColorimetricReference"],[50885,"SRawType"],[50898,"PanasonicTitle"],[50899,"PanasonicTitle2"],[50931,"CameraCalibrationSig"],[50932,"ProfileCalibrationSig"],[50933,"ProfileIFD"],[50934,"AsShotProfileName"],[50936,"ProfileName"],[50937,"ProfileHueSatMapDims"],[50938,"ProfileHueSatMapData1"],[50939,"ProfileHueSatMapData2"],[50940,"ProfileToneCurve"],[50941,"ProfileEmbedPolicy"],[50942,"ProfileCopyright"],[50964,"ForwardMatrix1"],[50965,"ForwardMatrix2"],[50966,"PreviewApplicationName"],[50967,"PreviewApplicationVersion"],[50968,"PreviewSettingsName"],[50969,"PreviewSettingsDigest"],[50970,"PreviewColorSpace"],[50971,"PreviewDateTime"],[50972,"RawImageDigest"],[50973,"OriginalRawFileDigest"],[50981,"ProfileLookTableDims"],[50982,"ProfileLookTableData"],[51043,"TimeCodes"],[51044,"FrameRate"],[51058,"TStop"],[51081,"ReelName"],[51089,"OriginalDefaultFinalSize"],[51090,"OriginalBestQualitySize"],[51091,"OriginalDefaultCropSize"],[51105,"CameraLabel"],[51107,"ProfileHueSatMapEncoding"],[51108,"ProfileLookTableEncoding"],[51109,"BaselineExposureOffset"],[51110,"DefaultBlackRender"],[51111,"NewRawImageDigest"],[51112,"RawToPreviewGain"]]);let Kr=[[273,"StripOffsets"],[279,"StripByteCounts"],[288,"FreeOffsets"],[289,"FreeByteCounts"],[291,"GrayResponseCurve"],[292,"T4Options"],[293,"T6Options"],[300,"ColorResponseUnit"],[320,"ColorMap"],[324,"TileOffsets"],[325,"TileByteCounts"],[326,"BadFaxLines"],[327,"CleanFaxData"],[328,"ConsecutiveBadFaxLines"],[330,"SubIFD"],[333,"InkNames"],[334,"NumberofInks"],[336,"DotRange"],[338,"ExtraSamples"],[339,"SampleFormat"],[340,"SMinSampleValue"],[341,"SMaxSampleValue"],[342,"TransferRange"],[343,"ClipPath"],[344,"XClipPathUnits"],[345,"YClipPathUnits"],[346,"Indexed"],[347,"JPEGTables"],[351,"OPIProxy"],[400,"GlobalParametersIFD"],[401,"ProfileType"],[402,"FaxProfile"],[403,"CodingMethods"],[404,"VersionYear"],[405,"ModeNumber"],[433,"Decode"],[434,"DefaultImageColor"],[435,"T82Options"],[437,"JPEGTables"],[512,"JPEGProc"],[515,"JPEGRestartInterval"],[517,"JPEGLosslessPredictors"],[518,"JPEGPointTransforms"],[519,"JPEGQTables"],[520,"JPEGDCTables"],[521,"JPEGACTables"],[559,"StripRowCounts"],[999,"USPTOMiscellaneous"],[18247,"XP_DIP_XML"],[18248,"StitchInfo"],[28672,"SonyRawFileType"],[28688,"SonyToneCurve"],[28721,"VignettingCorrection"],[28722,"VignettingCorrParams"],[28724,"ChromaticAberrationCorrection"],[28725,"ChromaticAberrationCorrParams"],[28726,"DistortionCorrection"],[28727,"DistortionCorrParams"],[29895,"SonyCropTopLeft"],[29896,"SonyCropSize"],[32781,"ImageID"],[32931,"WangTag1"],[32932,"WangAnnotation"],[32933,"WangTag3"],[32934,"WangTag4"],[32953,"ImageReferencePoints"],[32954,"RegionXformTackPoint"],[32955,"WarpQuadrilateral"],[32956,"AffineTransformMat"],[32995,"Matteing"],[32996,"DataType"],[32997,"ImageDepth"],[32998,"TileDepth"],[33300,"ImageFullWidth"],[33301,"ImageFullHeight"],[33302,"TextureFormat"],[33303,"WrapModes"],[33304,"FovCot"],[33305,"MatrixWorldToScreen"],[33306,"MatrixWorldToCamera"],[33405,"Model2"],[33421,"CFARepeatPatternDim"],[33422,"CFAPattern2"],[33423,"BatteryLevel"],[33424,"KodakIFD"],[33445,"MDFileTag"],[33446,"MDScalePixel"],[33447,"MDColorTable"],[33448,"MDLabName"],[33449,"MDSampleInfo"],[33450,"MDPrepDate"],[33451,"MDPrepTime"],[33452,"MDFileUnits"],[33589,"AdventScale"],[33590,"AdventRevision"],[33628,"UIC1Tag"],[33629,"UIC2Tag"],[33630,"UIC3Tag"],[33631,"UIC4Tag"],[33918,"IntergraphPacketData"],[33919,"IntergraphFlagRegisters"],[33921,"INGRReserved"],[34016,"Site"],[34017,"ColorSequence"],[34018,"IT8Header"],[34019,"RasterPadding"],[34020,"BitsPerRunLength"],[34021,"BitsPerExtendedRunLength"],[34022,"ColorTable"],[34023,"ImageColorIndicator"],[34024,"BackgroundColorIndicator"],[34025,"ImageColorValue"],[34026,"BackgroundColorValue"],[34027,"PixelIntensityRange"],[34028,"TransparencyIndicator"],[34029,"ColorCharacterization"],[34030,"HCUsage"],[34031,"TrapIndicator"],[34032,"CMYKEquivalent"],[34152,"AFCP_IPTC"],[34232,"PixelMagicJBIGOptions"],[34263,"JPLCartoIFD"],[34306,"WB_GRGBLevels"],[34310,"LeafData"],[34687,"TIFF_FXExtensions"],[34688,"MultiProfiles"],[34689,"SharedData"],[34690,"T88Options"],[34732,"ImageLayer"],[34750,"JBIGOptions"],[34856,"Opto-ElectricConvFactor"],[34857,"Interlace"],[34908,"FaxRecvParams"],[34909,"FaxSubAddress"],[34910,"FaxRecvTime"],[34929,"FedexEDR"],[34954,"LeafSubIFD"],[37387,"FlashEnergy"],[37388,"SpatialFrequencyResponse"],[37389,"Noise"],[37390,"FocalPlaneXResolution"],[37391,"FocalPlaneYResolution"],[37392,"FocalPlaneResolutionUnit"],[37397,"ExposureIndex"],[37398,"TIFF-EPStandardID"],[37399,"SensingMethod"],[37434,"CIP3DataFile"],[37435,"CIP3Sheet"],[37436,"CIP3Side"],[37439,"StoNits"],[37679,"MSDocumentText"],[37680,"MSPropertySetStorage"],[37681,"MSDocumentTextPosition"],[37724,"ImageSourceData"],[40965,"InteropIFD"],[40976,"SamsungRawPointersOffset"],[40977,"SamsungRawPointersLength"],[41217,"SamsungRawByteOrder"],[41218,"SamsungRawUnknown"],[41484,"SpatialFrequencyResponse"],[41485,"Noise"],[41489,"ImageNumber"],[41490,"SecurityClassification"],[41491,"ImageHistory"],[41494,"TIFF-EPStandardID"],[41995,"DeviceSettingDescription"],[42112,"GDALMetadata"],[42113,"GDALNoData"],[44992,"ExpandSoftware"],[44993,"ExpandLens"],[44994,"ExpandFilm"],[44995,"ExpandFilterLens"],[44996,"ExpandScanner"],[44997,"ExpandFlashLamp"],[46275,"HasselbladRawImage"],[48129,"PixelFormat"],[48130,"Transformation"],[48131,"Uncompressed"],[48132,"ImageType"],[48256,"ImageWidth"],[48257,"ImageHeight"],[48258,"WidthResolution"],[48259,"HeightResolution"],[48320,"ImageOffset"],[48321,"ImageByteCount"],[48322,"AlphaOffset"],[48323,"AlphaByteCount"],[48324,"ImageDataDiscard"],[48325,"AlphaDataDiscard"],[50215,"OceScanjobDesc"],[50216,"OceApplicationSelector"],[50217,"OceIDNumber"],[50218,"OceImageLogic"],[50255,"Annotations"],[50459,"HasselbladExif"],[50547,"OriginalFileName"],[50560,"USPTOOriginalContentType"],[50656,"CR2CFAPattern"],[50710,"CFAPlaneColor"],[50711,"CFALayout"],[50712,"LinearizationTable"],[50713,"BlackLevelRepeatDim"],[50714,"BlackLevel"],[50715,"BlackLevelDeltaH"],[50716,"BlackLevelDeltaV"],[50717,"WhiteLevel"],[50718,"DefaultScale"],[50719,"DefaultCropOrigin"],[50720,"DefaultCropSize"],[50733,"BayerGreenSplit"],[50737,"ChromaBlurRadius"],[50738,"AntiAliasStrength"],[50752,"RawImageSegmentation"],[50780,"BestQualityScale"],[50784,"AliasLayerMetadata"],[50829,"ActiveArea"],[50830,"MaskedAreas"],[50935,"NoiseReductionApplied"],[50974,"SubTileBlockSize"],[50975,"RowInterleaveFactor"],[51008,"OpcodeList1"],[51009,"OpcodeList2"],[51022,"OpcodeList3"],[51041,"NoiseProfile"],[51114,"CacheVersion"],[51125,"DefaultUserCrop"],[51157,"NikonNEFInfo"],[65024,"KdcIFD"]];on(kt,"ifd0",Kr),on(kt,"exif",Kr),wt($t,"gps",[[23,{M:"Magnetic North",T:"True North"}],[25,{K:"Kilometers",M:"Miles",N:"Nautical Miles"}]]);class Gi extends Kt{static canHandle(t,e){return t.getUint8(e+1)===224&&t.getUint32(e+4)===1246120262&&t.getUint8(e+8)===0}parse(){return this.parseTags(),this.translate(),this.output}parseTags(){this.raw=new Map([[0,this.chunk.getUint16(0)],[2,this.chunk.getUint8(2)],[3,this.chunk.getUint16(3)],[5,this.chunk.getUint16(5)],[7,this.chunk.getUint8(7)],[8,this.chunk.getUint8(8)]])}}tt(Gi,"type","jfif"),tt(Gi,"headerLength",9),Mt.set("jfif",Gi),wt(kt,"jfif",[[0,"JFIFVersion"],[2,"ResolutionUnit"],[3,"XResolution"],[5,"YResolution"],[7,"ThumbnailWidth"],[8,"ThumbnailHeight"]]);class $r extends Kt{parse(){return this.parseTags(),this.translate(),this.output}parseTags(){this.raw=new Map([[0,this.chunk.getUint32(0)],[4,this.chunk.getUint32(4)],[8,this.chunk.getUint8(8)],[9,this.chunk.getUint8(9)],[10,this.chunk.getUint8(10)],[11,this.chunk.getUint8(11)],[12,this.chunk.getUint8(12)],...Array.from(this.raw)])}}tt($r,"type","ihdr"),Mt.set("ihdr",$r),wt(kt,"ihdr",[[0,"ImageWidth"],[4,"ImageHeight"],[8,"BitDepth"],[9,"ColorType"],[10,"Compression"],[11,"Filter"],[12,"Interlace"]]),wt($t,"ihdr",[[9,{0:"Grayscale",2:"RGB",3:"Palette",4:"Grayscale with Alpha",6:"RGB with Alpha",DEFAULT:"Unknown"}],[10,{0:"Deflate/Inflate",DEFAULT:"Unknown"}],[11,{0:"Adaptive",DEFAULT:"Unknown"}],[12,{0:"Noninterlaced",1:"Adam7 Interlace",DEFAULT:"Unknown"}]]);class ti extends Kt{static canHandle(t,e){return t.getUint8(e+1)===226&&t.getUint32(e+4)===1229144927}static findPosition(t,e){let i=super.findPosition(t,e);return i.chunkNumber=t.getUint8(e+16),i.chunkCount=t.getUint8(e+17),i.multiSegment=i.chunkCount>1,i}static handleMultiSegments(t){return function(e){let i=function(r){let s=r[0].constructor,o=0;for(let u of r)o+=u.length;let a=new s(o),l=0;for(let u of r)a.set(u,l),l+=u.length;return a}(e.map(r=>r.chunk.toUint8()));return new zt(i)}(t)}parse(){return this.raw=new Map,this.parseHeader(),this.parseTags(),this.translate(),this.output}parseHeader(){let{raw:t}=this;this.chunk.byteLength<84&&Pt("ICC header is too short");for(let[e,i]of Object.entries(Za)){e=parseInt(e,10);let r=i(this.chunk,e);r!=="\0\0\0\0"&&t.set(e,r)}}parseTags(){let t,e,i,r,s,{raw:o}=this,a=this.chunk.getUint32(128),l=132,u=this.chunk.byteLength;for(;a--;){if(t=this.chunk.getString(l,4),e=this.chunk.getUint32(l+4),i=this.chunk.getUint32(l+8),r=this.chunk.getString(e,4),e+i>u)return void console.warn("reached the end of the first ICC chunk. Enable options.tiff.multiSegment to read all ICC segments.");s=this.parseTag(r,e,i),s!==void 0&&s!=="\0\0\0\0"&&o.set(t,s),l+=12}}parseTag(t,e,i){switch(t){case"desc":return this.parseDesc(e);case"mluc":return this.parseMluc(e);case"text":return this.parseText(e,i);case"sig ":return this.parseSig(e)}if(!(e+i>this.chunk.byteLength))return this.chunk.getUint8Array(e,i)}parseDesc(t){let e=this.chunk.getUint32(t+8)-1;return Je(this.chunk.getString(t+12,e))}parseText(t,e){return Je(this.chunk.getString(t+8,e-8))}parseSig(t){return Je(this.chunk.getString(t+8,4))}parseMluc(t){let{chunk:e}=this,i=e.getUint32(t+8),r=e.getUint32(t+12),s=t+16,o=[];for(let a=0;a<i;a++){let l=e.getString(s+0,2),u=e.getString(s+2,2),h=e.getUint32(s+4),c=e.getUint32(s+8)+t,m=Je(e.getUnicodeString(c,h));o.push({lang:l,country:u,text:m}),s+=r}return i===1?o[0].text:o}translateValue(t,e){return typeof t=="string"?e[t]||e[t.toLowerCase()]||t:e[t]||t}}tt(ti,"type","icc"),tt(ti,"multiSegment",!0),tt(ti,"headerLength",18);const Za={4:ye,8:function(n,t){return[n.getUint8(t),n.getUint8(t+1)>>4,n.getUint8(t+1)%16].map(e=>e.toString(10)).join(".")},12:ye,16:ye,20:ye,24:function(n,t){const e=n.getUint16(t),i=n.getUint16(t+2)-1,r=n.getUint16(t+4),s=n.getUint16(t+6),o=n.getUint16(t+8),a=n.getUint16(t+10);return new Date(Date.UTC(e,i,r,s,o,a))},36:ye,40:ye,48:ye,52:ye,64:(n,t)=>n.getUint32(t),80:ye};function ye(n,t){return Je(n.getString(t,4))}Mt.set("icc",ti),wt(kt,"icc",[[4,"ProfileCMMType"],[8,"ProfileVersion"],[12,"ProfileClass"],[16,"ColorSpaceData"],[20,"ProfileConnectionSpace"],[24,"ProfileDateTime"],[36,"ProfileFileSignature"],[40,"PrimaryPlatform"],[44,"CMMFlags"],[48,"DeviceManufacturer"],[52,"DeviceModel"],[56,"DeviceAttributes"],[64,"RenderingIntent"],[68,"ConnectionSpaceIlluminant"],[80,"ProfileCreator"],[84,"ProfileID"],["Header","ProfileHeader"],["MS00","WCSProfiles"],["bTRC","BlueTRC"],["bXYZ","BlueMatrixColumn"],["bfd","UCRBG"],["bkpt","MediaBlackPoint"],["calt","CalibrationDateTime"],["chad","ChromaticAdaptation"],["chrm","Chromaticity"],["ciis","ColorimetricIntentImageState"],["clot","ColorantTableOut"],["clro","ColorantOrder"],["clrt","ColorantTable"],["cprt","ProfileCopyright"],["crdi","CRDInfo"],["desc","ProfileDescription"],["devs","DeviceSettings"],["dmdd","DeviceModelDesc"],["dmnd","DeviceMfgDesc"],["dscm","ProfileDescriptionML"],["fpce","FocalPlaneColorimetryEstimates"],["gTRC","GreenTRC"],["gXYZ","GreenMatrixColumn"],["gamt","Gamut"],["kTRC","GrayTRC"],["lumi","Luminance"],["meas","Measurement"],["meta","Metadata"],["mmod","MakeAndModel"],["ncl2","NamedColor2"],["ncol","NamedColor"],["ndin","NativeDisplayInfo"],["pre0","Preview0"],["pre1","Preview1"],["pre2","Preview2"],["ps2i","PS2RenderingIntent"],["ps2s","PostScript2CSA"],["psd0","PostScript2CRD0"],["psd1","PostScript2CRD1"],["psd2","PostScript2CRD2"],["psd3","PostScript2CRD3"],["pseq","ProfileSequenceDesc"],["psid","ProfileSequenceIdentifier"],["psvm","PS2CRDVMSize"],["rTRC","RedTRC"],["rXYZ","RedMatrixColumn"],["resp","OutputResponse"],["rhoc","ReflectionHardcopyOrigColorimetry"],["rig0","PerceptualRenderingIntentGamut"],["rig2","SaturationRenderingIntentGamut"],["rpoc","ReflectionPrintOutputColorimetry"],["sape","SceneAppearanceEstimates"],["scoe","SceneColorimetryEstimates"],["scrd","ScreeningDesc"],["scrn","Screening"],["targ","CharTarget"],["tech","Technology"],["vcgt","VideoCardGamma"],["view","ViewingConditions"],["vued","ViewingCondDesc"],["wtpt","MediaWhitePoint"]]);const jn={"4d2p":"Erdt Systems",AAMA:"Aamazing Technologies",ACER:"Acer",ACLT:"Acolyte Color Research",ACTI:"Actix Sytems",ADAR:"Adara Technology",ADBE:"Adobe",ADI:"ADI Systems",AGFA:"Agfa Graphics",ALMD:"Alps Electric",ALPS:"Alps Electric",ALWN:"Alwan Color Expertise",AMTI:"Amiable Technologies",AOC:"AOC International",APAG:"Apago",APPL:"Apple Computer",AST:"AST","AT&T":"AT&T",BAEL:"BARBIERI electronic",BRCO:"Barco NV",BRKP:"Breakpoint",BROT:"Brother",BULL:"Bull",BUS:"Bus Computer Systems","C-IT":"C-Itoh",CAMR:"Intel",CANO:"Canon",CARR:"Carroll Touch",CASI:"Casio",CBUS:"Colorbus PL",CEL:"Crossfield",CELx:"Crossfield",CGS:"CGS Publishing Technologies International",CHM:"Rochester Robotics",CIGL:"Colour Imaging Group, London",CITI:"Citizen",CL00:"Candela",CLIQ:"Color IQ",CMCO:"Chromaco",CMiX:"CHROMiX",COLO:"Colorgraphic Communications",COMP:"Compaq",COMp:"Compeq/Focus Technology",CONR:"Conrac Display Products",CORD:"Cordata Technologies",CPQ:"Compaq",CPRO:"ColorPro",CRN:"Cornerstone",CTX:"CTX International",CVIS:"ColorVision",CWC:"Fujitsu Laboratories",DARI:"Darius Technology",DATA:"Dataproducts",DCP:"Dry Creek Photo",DCRC:"Digital Contents Resource Center, Chung-Ang University",DELL:"Dell Computer",DIC:"Dainippon Ink and Chemicals",DICO:"Diconix",DIGI:"Digital","DL&C":"Digital Light & Color",DPLG:"Doppelganger",DS:"Dainippon Screen",DSOL:"DOOSOL",DUPN:"DuPont",EPSO:"Epson",ESKO:"Esko-Graphics",ETRI:"Electronics and Telecommunications Research Institute",EVER:"Everex Systems",EXAC:"ExactCODE",Eizo:"Eizo",FALC:"Falco Data Products",FF:"Fuji Photo Film",FFEI:"FujiFilm Electronic Imaging",FNRD:"Fnord Software",FORA:"Fora",FORE:"Forefront Technology",FP:"Fujitsu",FPA:"WayTech Development",FUJI:"Fujitsu",FX:"Fuji Xerox",GCC:"GCC Technologies",GGSL:"Global Graphics Software",GMB:"Gretagmacbeth",GMG:"GMG",GOLD:"GoldStar Technology",GOOG:"Google",GPRT:"Giantprint",GTMB:"Gretagmacbeth",GVC:"WayTech Development",GW2K:"Sony",HCI:"HCI",HDM:"Heidelberger Druckmaschinen",HERM:"Hermes",HITA:"Hitachi America",HP:"Hewlett-Packard",HTC:"Hitachi",HiTi:"HiTi Digital",IBM:"IBM",IDNT:"Scitex",IEC:"Hewlett-Packard",IIYA:"Iiyama North America",IKEG:"Ikegami Electronics",IMAG:"Image Systems",IMI:"Ingram Micro",INTC:"Intel",INTL:"N/A (INTL)",INTR:"Intra Electronics",IOCO:"Iocomm International Technology",IPS:"InfoPrint Solutions Company",IRIS:"Scitex",ISL:"Ichikawa Soft Laboratory",ITNL:"N/A (ITNL)",IVM:"IVM",IWAT:"Iwatsu Electric",Idnt:"Scitex",Inca:"Inca Digital Printers",Iris:"Scitex",JPEG:"Joint Photographic Experts Group",JSFT:"Jetsoft Development",JVC:"JVC Information Products",KART:"Scitex",KFC:"KFC Computek Components",KLH:"KLH Computers",KMHD:"Konica Minolta",KNCA:"Konica",KODA:"Kodak",KYOC:"Kyocera",Kart:"Scitex",LCAG:"Leica",LCCD:"Leeds Colour",LDAK:"Left Dakota",LEAD:"Leading Technology",LEXM:"Lexmark International",LINK:"Link Computer",LINO:"Linotronic",LITE:"Lite-On",Leaf:"Leaf",Lino:"Linotronic",MAGC:"Mag Computronic",MAGI:"MAG Innovision",MANN:"Mannesmann",MICN:"Micron Technology",MICR:"Microtek",MICV:"Microvitec",MINO:"Minolta",MITS:"Mitsubishi Electronics America",MITs:"Mitsuba",MNLT:"Minolta",MODG:"Modgraph",MONI:"Monitronix",MONS:"Monaco Systems",MORS:"Morse Technology",MOTI:"Motive Systems",MSFT:"Microsoft",MUTO:"MUTOH INDUSTRIES",Mits:"Mitsubishi Electric",NANA:"NANAO",NEC:"NEC",NEXP:"NexPress Solutions",NISS:"Nissei Sangyo America",NKON:"Nikon",NONE:"none",OCE:"Oce Technologies",OCEC:"OceColor",OKI:"Oki",OKID:"Okidata",OKIP:"Okidata",OLIV:"Olivetti",OLYM:"Olympus",ONYX:"Onyx Graphics",OPTI:"Optiquest",PACK:"Packard Bell",PANA:"Matsushita Electric Industrial",PANT:"Pantone",PBN:"Packard Bell",PFU:"PFU",PHIL:"Philips Consumer Electronics",PNTX:"HOYA",POne:"Phase One A/S",PREM:"Premier Computer Innovations",PRIN:"Princeton Graphic Systems",PRIP:"Princeton Publishing Labs",QLUX:"Hong Kong",QMS:"QMS",QPCD:"QPcard AB",QUAD:"QuadLaser",QUME:"Qume",RADI:"Radius",RDDx:"Integrated Color Solutions",RDG:"Roland DG",REDM:"REDMS Group",RELI:"Relisys",RGMS:"Rolf Gierling Multitools",RICO:"Ricoh",RNLD:"Edmund Ronald",ROYA:"Royal",RPC:"Ricoh Printing Systems",RTL:"Royal Information Electronics",SAMP:"Sampo",SAMS:"Samsung",SANT:"Jaime Santana Pomares",SCIT:"Scitex",SCRN:"Dainippon Screen",SDP:"Scitex",SEC:"Samsung",SEIK:"Seiko Instruments",SEIk:"Seikosha",SGUY:"ScanGuy.com",SHAR:"Sharp Laboratories",SICC:"International Color Consortium",SONY:"Sony",SPCL:"SpectraCal",STAR:"Star",STC:"Sampo Technology",Scit:"Scitex",Sdp:"Scitex",Sony:"Sony",TALO:"Talon Technology",TAND:"Tandy",TATU:"Tatung",TAXA:"TAXAN America",TDS:"Tokyo Denshi Sekei",TECO:"TECO Information Systems",TEGR:"Tegra",TEKT:"Tektronix",TI:"Texas Instruments",TMKR:"TypeMaker",TOSB:"Toshiba",TOSH:"Toshiba",TOTK:"TOTOKU ELECTRIC",TRIU:"Triumph",TSBT:"Toshiba",TTX:"TTX Computer Products",TVM:"TVM Professional Monitor",TW:"TW Casper",ULSX:"Ulead Systems",UNIS:"Unisys",UTZF:"Utz Fehlau & Sohn",VARI:"Varityper",VIEW:"Viewsonic",VISL:"Visual communication",VIVO:"Vivo Mobile Communication",WANG:"Wang",WLBR:"Wilbur Imaging",WTG2:"Ware To Go",WYSE:"WYSE Technology",XERX:"Xerox",XRIT:"X-Rite",ZRAN:"Zoran",Zebr:"Zebra Technologies",appl:"Apple Computer",bICC:"basICColor",berg:"bergdesign",ceyd:"Integrated Color Solutions",clsp:"MacDermid ColorSpan",ds:"Dainippon Screen",dupn:"DuPont",ffei:"FujiFilm Electronic Imaging",flux:"FluxData",iris:"Scitex",kart:"Scitex",lcms:"Little CMS",lino:"Linotronic",none:"none",ob4d:"Erdt Systems",obic:"Medigraph",quby:"Qubyx Sarl",scit:"Scitex",scrn:"Dainippon Screen",sdp:"Scitex",siwi:"SIWI GRAFIKA",yxym:"YxyMaster"},Jr={scnr:"Scanner",mntr:"Monitor",prtr:"Printer",link:"Device Link",abst:"Abstract",spac:"Color Space Conversion Profile",nmcl:"Named Color",cenc:"ColorEncodingSpace profile",mid:"MultiplexIdentification profile",mlnk:"MultiplexLink profile",mvis:"MultiplexVisualization profile",nkpf:"Nikon Input Device Profile (NON-STANDARD!)"};wt($t,"icc",[[4,jn],[12,Jr],[40,Object.assign({},jn,Jr)],[48,jn],[80,jn],[64,{0:"Perceptual",1:"Relative Colorimetric",2:"Saturation",3:"Absolute Colorimetric"}],["tech",{amd:"Active Matrix Display",crt:"Cathode Ray Tube Display",kpcd:"Photo CD",pmd:"Passive Matrix Display",dcam:"Digital Camera",dcpj:"Digital Cinema Projector",dmpc:"Digital Motion Picture Camera",dsub:"Dye Sublimation Printer",epho:"Electrophotographic Printer",esta:"Electrostatic Printer",flex:"Flexography",fprn:"Film Writer",fscn:"Film Scanner",grav:"Gravure",ijet:"Ink Jet Printer",imgs:"Photo Image Setter",mpfr:"Motion Picture Film Recorder",mpfs:"Motion Picture Film Scanner",offs:"Offset Lithography",pjtv:"Projection Television",rpho:"Photographic Paper Printer",rscn:"Reflective Scanner",silk:"Silkscreen",twax:"Thermal Wax Printer",vidc:"Video Camera",vidm:"Video Monitor"}]]);class qn extends Kt{static canHandle(t,e,i){return t.getUint8(e+1)===237&&t.getString(e+4,9)==="Photoshop"&&this.containsIptc8bim(t,e,i)!==void 0}static headerLength(t,e,i){let r,s=this.containsIptc8bim(t,e,i);if(s!==void 0)return r=t.getUint8(e+s+7),r%2!=0&&(r+=1),r===0&&(r=4),s+8+r}static containsIptc8bim(t,e,i){for(let r=0;r<i;r++)if(this.isIptcSegmentHead(t,e+r))return r}static isIptcSegmentHead(t,e){return t.getUint8(e)===56&&t.getUint32(e)===943868237&&t.getUint16(e+4)===1028}parse(){let{raw:t}=this,e=this.chunk.byteLength-1,i=!1;for(let r=0;r<e;r++)if(this.chunk.getUint8(r)===28&&this.chunk.getUint8(r+1)===2){i=!0;let s=this.chunk.getUint16(r+3),o=this.chunk.getUint8(r+2),a=this.chunk.getLatin1String(r+5,s);t.set(o,this.pluralizeValue(t.get(o),a)),r+=4+s}else if(i)break;return this.translate(),this.output}pluralizeValue(t,e){return t!==void 0?t instanceof Array?(t.push(e),t):[t,e]:e}}tt(qn,"type","iptc"),tt(qn,"translateValues",!1),tt(qn,"reviveValues",!1),Mt.set("iptc",qn),wt(kt,"iptc",[[0,"ApplicationRecordVersion"],[3,"ObjectTypeReference"],[4,"ObjectAttributeReference"],[5,"ObjectName"],[7,"EditStatus"],[8,"EditorialUpdate"],[10,"Urgency"],[12,"SubjectReference"],[15,"Category"],[20,"SupplementalCategories"],[22,"FixtureIdentifier"],[25,"Keywords"],[26,"ContentLocationCode"],[27,"ContentLocationName"],[30,"ReleaseDate"],[35,"ReleaseTime"],[37,"ExpirationDate"],[38,"ExpirationTime"],[40,"SpecialInstructions"],[42,"ActionAdvised"],[45,"ReferenceService"],[47,"ReferenceDate"],[50,"ReferenceNumber"],[55,"DateCreated"],[60,"TimeCreated"],[62,"DigitalCreationDate"],[63,"DigitalCreationTime"],[65,"OriginatingProgram"],[70,"ProgramVersion"],[75,"ObjectCycle"],[80,"Byline"],[85,"BylineTitle"],[90,"City"],[92,"Sublocation"],[95,"State"],[100,"CountryCode"],[101,"Country"],[103,"OriginalTransmissionReference"],[105,"Headline"],[110,"Credit"],[115,"Source"],[116,"CopyrightNotice"],[118,"Contact"],[120,"Caption"],[121,"LocalCaption"],[122,"Writer"],[125,"RasterizedCaption"],[130,"ImageType"],[131,"ImageOrientation"],[135,"LanguageIdentifier"],[150,"AudioType"],[151,"AudioSamplingRate"],[152,"AudioSamplingResolution"],[153,"AudioDuration"],[154,"AudioOutcue"],[184,"JobID"],[185,"MasterDocumentID"],[186,"ShortDocumentID"],[187,"UniqueDocumentID"],[188,"OwnerID"],[200,"ObjectPreviewFileFormat"],[201,"ObjectPreviewFileVersion"],[202,"ObjectPreviewData"],[221,"Prefs"],[225,"ClassifyState"],[228,"SimilarityIndex"],[230,"DocumentNotes"],[231,"DocumentHistory"],[232,"ExifCameraInfo"],[255,"CatalogSets"]]),wt($t,"iptc",[[10,{0:"0 (reserved)",1:"1 (most urgent)",2:"2",3:"3",4:"4",5:"5 (normal urgency)",6:"6",7:"7",8:"8 (least urgent)",9:"9 (user-defined priority)"}],[75,{a:"Morning",b:"Both Morning and Evening",p:"Evening"}],[131,{L:"Landscape",P:"Portrait",S:"Square"}]]);function to(n,t,e){const i=e.type===1||e.type===2||e.type===7?1:e.type===3?2:e.type===4||e.type===9||e.type===11||e.type===13?4:e.type===5||e.type===10||e.type===12?8:0;if(i===0)throw new Error(`Unsupported TIFF field type: ${e.type}`);const s=e.count*i<=4?e.valueFieldOffset:n.getUint32(e.valueFieldOffset,t),o=[];for(let a=0;a<e.count;a++){const l=s+a*i;if(l<0||l+i>n.byteLength)throw new Error("Invalid TIFF field offset");if(e.type===1)o.push(n.getUint8(l));else if(e.type===2||e.type===7)o.push(n.getUint8(l));else if(e.type===3)o.push(n.getUint16(l,t));else if(e.type===5){const u=n.getUint32(l+4,t);o.push(u?n.getUint32(l,t)/u:0)}else if(e.type===13)o.push(n.getUint32(l,t));else if(e.type===9)o.push(n.getInt32(l,t));else if(e.type===10){const u=n.getInt32(l+4,t);o.push(u?n.getInt32(l,t)/u:0)}else e.type===11?o.push(n.getFloat32(l,t)):e.type===12?o.push(n.getFloat64(l,t)):o.push(n.getUint32(l,t))}return o}function Zr(n,t,e){if(!Number.isFinite(e)||e<=0||e+2>n.byteLength)throw new Error("Invalid TIFF IFD offset");const i=n.getUint16(e,t);if(e+2+i*12+4>n.byteLength)throw new Error("Corrupt TIFF IFD");const s=new Map;for(let o=0;o<i;o++){const a=e+2+o*12,l=n.getUint16(a,t),u=n.getUint16(a+2,t),h=n.getUint32(a+4,t);s.set(l,{type:u,count:h,valueFieldOffset:a+8})}return s}function Ot(n,t,e,i){const r=e.get(i);return r?to(n,t,r):[]}const Ts=(...n)=>{for(const t of n){if(Array.isArray(t)||ArrayBuffer.isView(t)){const i=Array.from(t).map(Number).filter(r=>Number.isFinite(r)&&r>0);if(i.length)return Math.max(...i);continue}const e=Number(t);if(Number.isFinite(e)&&e>0)return e}return null},eo=n=>{const e=(Array.isArray(n)||ArrayBuffer.isView(n)?Array.from(n).map(Number):[Number(n)]).map(i=>Number.isFinite(i)?i:0);return e.length>=3?[e[0]||0,e[1]||0,e[1]||0,e[2]||0]:e.length===1?[e[0]||0,e[0]||0,e[0]||0,e[0]||0]:[0,0,0,0]};function no(n,t={}){var v,k,A,S,F,T;if(n.byteLength<8)return null;const e=new DataView(n.buffer,n.byteOffset,n.byteLength),i=e.getUint16(0,!1),r=i===18761;if(!r&&i!==19789||e.getUint16(2,r)!==42)return null;const s=e.getUint32(4,r);let o;try{o=Zr(e,r,s)}catch{return null}const a=new Set([s]);for(const R of Ot(e,r,o,330))R>0&&a.add(R);let l=null;for(const R of a)try{const E=R===s?o:Zr(e,r,R),O=Ot(e,r,E,256)[0]||0,U=Ot(e,r,E,257)[0]||0,I=Ot(e,r,E,258),D=Ot(e,r,E,277)[0]||I.length||1,B=I.length===1?new Array(D).fill(I[0]):I,V=Ot(e,r,E,259)[0]||1,Q=Ot(e,r,E,262)[0]||0,z=Ot(e,r,E,284)[0]||1,Y=Ot(e,r,E,339),W=Y.length===0?new Array(D).fill(1):Y.length===1?new Array(D).fill(Y[0]):Y,K=Ot(e,r,E,273),it=Ot(e,r,E,279),J=Ot(e,r,E,278)[0]||U,ot=W.slice(0,D).every(L=>L===0||L===1);if(!(Q===34892&&V===1&&z===1&&O>0&&U>0&&D>=3&&B.length>=D&&B.slice(0,D).every(L=>L===16)&&ot&&K.length>0&&K.length===it.length))continue;(!l||O*U>l.width*l.height)&&(l={entries:E,width:O,height:U,samplesPerPixel:D,rowsPerStrip:J,stripOffsets:K,stripByteCounts:it,bitsPerSample:B,photometric:Q})}catch{continue}if(!l)return null;const{entries:u,width:h,height:c,samplesPerPixel:m,rowsPerStrip:d,stripOffsets:p,stripByteCounts:f,bitsPerSample:y,photometric:x}=l,g=h*c;if(!Number.isSafeInteger(g)||g<=0)return null;const b=new Uint16Array(g*3);for(let R=0;R<p.length;R++){const E=R*d;if(E>=c)break;const O=Math.min(d,c-E),U=O*h*m*2,I=p[R],D=f[R];if(I<0||D<U||I+U>n.byteLength)throw new Error("Invalid LinearRaw DNG strip bounds");const B=E*h*3;if(m===3&&r&&!(n.byteOffset+I&1)){b.set(new Uint16Array(n.buffer,n.byteOffset+I,O*h*3),B);continue}const V=new DataView(n.buffer,n.byteOffset+I,U);let Q=0,z=B;for(let Y=0;Y<O*h;Y++)b[z++]=V.getUint16(Q,r),b[z++]=V.getUint16(Q+2,r),b[z++]=V.getUint16(Q+4,r),Q+=m*2}const _=Ot(e,r,u,50714),M=((k=(v=t==null?void 0:t.color_data)==null?void 0:v.dng_levels)==null?void 0:k.dng_cblack)||((A=t==null?void 0:t.color_data)==null?void 0:A.cblack_rawpy_style)||(t==null?void 0:t.black_level_per_channel)||(t==null?void 0:t.cblack),w=eo(_.length?_:M),P=Ot(e,r,u,50717),C=Ts(P,(F=(S=t==null?void 0:t.color_data)==null?void 0:S.dng_levels)==null?void 0:F.dng_whitelevel,(T=t==null?void 0:t.color_data)==null?void 0:T.maximum,t==null?void 0:t.white_level)||65535;return{data:b,width:h,height:c,bayerPattern:"",blackLevels:w,whiteLevel:C,metadata:{...t,format:"DNG_LINEAR_RAW_RGB",description:"Uncompressed DNG LinearRaw RGB",linearRawDngDecoder:!0,bitsPerSample:y.slice(0,m),samplesPerPixel:m,photometric:x},isThreePlane:!0,threePlaneTransfer:"linear"}}async function io(n){const t=no(n);if(!t)return null;let e={};const i=await tr(),r=new i;try{await r.open(n,{}),e=await r.metadata(!0)}catch(s){console.warn("LinearRaw DNG metadata enrichment failed",s)}finally{try{r.delete?r.delete():r.close()}catch{}}try{const s=await As.parse(n.buffer);s&&(e={...e,...s})}catch(s){console.warn("exifr parsing failed for LinearRaw DNG",s)}return t.metadata={...e,...t.metadata},t}const ro=async n=>{var r,s,o,a,l,u,h,c,m,d,p;const t=await io(n);if(t)return t;const e=await tr(),i=new e;try{if(await i.open(n,{}),typeof i.getRawImage!="function")throw new Error("WASM mismatch");let f={};try{f=await i.metadata(!0)}catch(F){console.warn("Metadata error before raw extraction",F)}const y=i.getRawImage(),x=y.data instanceof Uint16Array?y.data:new Uint16Array(y.data);let g={...f};try{const F=await As.parse(n.buffer);F&&(g={...g,...F})}catch(F){console.warn("exifr parsing failed for RAW buffer",F)}const b=f.filters||((r=f.idata)==null?void 0:r.filters)||0,_=f.colors||((s=f.idata)==null?void 0:s.colors)||0,M=b===0&&_===3,w=b===9;let P=[0,0,0,0],C=null,v;if(i.getBlackLevels)try{const F=i.getBlackLevels();C=F,v=Lr(F)||Lr((o=f.color_data)==null?void 0:o.black_level_model)||void 0;const T=Ke(v==null?void 0:v.siteBaseLevels);T&&(P=T)}catch(F){console.warn("getBlackLevels binding failed",F)}if(!v){const F=Number((C==null?void 0:C.black)??((a=f.color_data)==null?void 0:a.black)??0)||0,T=Ke((C==null?void 0:C.cblack)||((l=f.color_data)==null?void 0:l.cblack_rawpy_style)||f.black_level_per_channel||f.cblack||((u=f.color)==null?void 0:u.cblack));if(T){const R=String(y.bayerPattern||f.cfa_pattern||"RGGB").toUpperCase();let E=0;P=[0,1,2,3].map(O=>{const U=R[O],I=U==="R"?0:U==="B"?2:++E===1?1:3;return Math.max(0,F+T[I])})}else P=[F,F,F,F]}const k=Number.isFinite(Number(y.bits))&&Number(y.bits)>0?Math.pow(2,Number(y.bits))-1:null,A=Ts(f.white_level,(c=(h=f.color_data)==null?void 0:h.dng_levels)==null?void 0:c.dng_whitelevel,(m=f.color_data)==null?void 0:m.maximum,C==null?void 0:C.maximum,(d=f.color_data)==null?void 0:d.fmaximum,(p=f.color_data)==null?void 0:p.data_maximum,k)||16383,S={data:x,width:y.width,height:y.height,bayerPattern:y.bayerPattern||"",blackLevels:P,blackLevelModel:v,whiteLevel:A,metadata:g,isThreePlane:M,threePlaneTransfer:M?"linear":void 0,isXTrans:w};return ka(S,Fa(S,g)),S}finally{i.delete?i.delete():i.close()}};async function so(n){if(xa(n)){const i=await Sa(n);if(!i)throw new Error("Sony cRAW HQ decoder did not return image data.");return i.rawImageData}const e=new Uint8Array(n);return ro(e)}function ao(n,t,e){const i=Is(n,e),r=Math.floor(t.x),s=Math.floor(t.y),o=Math.floor(t.w),a=Math.floor(t.h),l=new Uint16Array(o*a);for(let u=0;u<a;u++){const h=s+u,c=u*o;for(let m=0;m<o;m++)l[c+m]=i(r+m,h)}return{data:l,width:o,height:a}}function oo(n,t,e){const i=lo(n,t,e);if(!i)return null;const r=n.width,s=n.height,o=new Uint16Array(r*s);for(let a=0;a<s;a++){const l=a*r;for(let u=0;u<r;u++)o[l+u]=i(u,a)}return{kind:"u16-mono",data:o,width:r,height:s}}function lo(n,t,e){return t.renderMode==="advanced-zero-dep"&&t.advancedZeroDep?Is(n,t,e):t.renderMode==="zero-dependency"?co(n,t,e):null}function Is(n,t,e){if(!t.advancedZeroDep)throw new Error("Unmixing settings not found in DisplaySettings.");const{bg:i,fg:r}=t.advancedZeroDep,s=Ns(e,t.advancedZeroDep.bl),{data:o,width:a,whiteLevel:l}=n,u=i.map((p,f)=>Math.max(0,p-s[f])),h=r.map((p,f)=>Math.max(0,p-s[f])),c=(u[1]+u[3])/2,m=(h[1]+h[3])/2,d=Math.pow(2,t.exposure);return(p,f)=>{if(p<0||f<0||p>=a||f>=n.height)return 0;const y=f%2,x=p%2;let g=0;!y&&!x?g=0:!y&&x?g=1:y&&!x?g=3:g=2;const b=o[f*a+p],_=s[g],M=u[g],w=h[g],P=Math.max(b-_,0),C=w-M||1e-9,v=(P-M)/C;let k;return v<0?k=P*(c/Math.max(M,1e-9)):v>1?k=P*(m/Math.max(w,1e-9)):k=(1-v)*c+v*m,k*=d,Math.max(0,Math.min(65535,Math.round(k)))}}function co(n,t,e){const{data:i,width:r,height:s,whiteLevel:o}=n,a=Ns(e,t.blackLevel||[0,0,0,0]),l=ho(n.bayerPattern),u=t.wbGains?t.wbGains[0]:1,h=t.wbGains?t.wbGains[1]:1,c=Math.pow(2,t.exposure||0);return(m,d)=>{if(m<0||d<0||m>=r||d>=s)return 0;const p=fo(m,d),f=po(l,m,d),y=i[d*r+m],x=a[p];let g=(y-x)/Math.max(1,o-x);return g=Math.max(0,Math.min(1,g)),g*=c,f==="R"?g*=u:f==="B"&&(g*=h),uo(g)}}function Ns(n,t){if(typeof n=="number"&&Number.isFinite(n)){const e=Math.max(0,n);return[e,e,e,e]}return Array.isArray(n)&&n.length===4?[Math.max(0,n[0]??0),Math.max(0,n[1]??0),Math.max(0,n[2]??0),Math.max(0,n[3]??0)]:[Math.max(0,t[0]??0),Math.max(0,t[1]??0),Math.max(0,t[2]??0),Math.max(0,t[3]??0)]}function uo(n){return Math.max(0,Math.min(65535,Math.round(Math.max(0,Math.min(1,n))*65535)))}function ho(n){const t=$e(n);if(!t)throw new Error("Cannot process RAW mosaic without a valid Bayer CFA pattern.");return t}function fo(n,t){return ms(n,t)}function po(n,t,e){return Aa(n,t,e)}const ut=(n,t=0)=>({real:n,imag:t}),mn=(n,t)=>({real:n.real+t.real,imag:n.imag+t.imag}),dn=(n,t)=>({real:n.real-t.real,imag:n.imag-t.imag}),Ct=(n,t)=>({real:n.real*t.real-n.imag*t.imag,imag:n.real*t.imag+n.imag*t.real}),Te=(n,t)=>{const e=t.real*t.real+t.imag*t.imag;return e===0?ut(0):{real:(n.real*t.real+n.imag*t.imag)/e,imag:(n.imag*t.real-n.real*t.imag)/e}},Qe=n=>Math.hypot(n.real,n.imag),Rs=n=>{const t=Qe(n);if(t===0)return ut(0);const e=Math.sqrt(t),i=Math.atan2(n.imag,n.real);return ut(e*Math.cos(i/2),e*Math.sin(i/2))};function mo(n,t){const e=n.length-1;if(e<0)return{p:ut(0),dp:ut(0),d2p:ut(0)};let i=ut(n[e].real,n[e].imag),r=ut(0),s=ut(0);for(let o=e-1;o>=0;o--)s=mn(Ct(r,ut(2)),Ct(t,s)),r=mn(i,Ct(t,r)),i=mn(ut(n[o].real,n[o].imag),Ct(t,i));return{p:i,dp:r,d2p:s}}function $i(n,t,e=80){const r=n.length-1;if(r<=0)return{root:t,iterations:0};if(r===1)return{root:Te(Ct(n[0],ut(-1)),n[1]),iterations:0};let s=ut(t.real,t.imag);for(let o=0;o<e;o++){const{p:a,dp:l,d2p:u}=mo(n,s);if(Qe(a)<1e-14)return{root:s,iterations:o};const h=Te(l,a),c=Ct(h,h),m=dn(c,Te(u,a)),d=ut(r),p=ut(r-1),f=dn(Ct(d,m),Ct(h,h)),y=Rs(Ct(p,f)),x=mn(h,y),g=dn(h,y),b=Qe(x)>Qe(g)?x:g;if(Qe(b)<1e-14)return{root:s,iterations:o};const _=Te(d,b),M=dn(s,_);if(Qe(_)<1e-14*Qe(M))return{root:M,iterations:o+1};s=M}return{root:s,iterations:e}}function go(n,t){const e=n.length-1;if(e<=0)return[ut(0)];if(e===1)return[ut(n[0].real,n[0].imag)];const i=new Array(e);i[e-1]=ut(n[e].real,n[e].imag);for(let r=e-2;r>=0;r--){const s=ut(n[r+1].real,n[r+1].imag),o=i[r+1];i[r]=mn(s,Ct(t,o))}return i}function yo(n){const t=n.length-1;if(t<=0)return[];if(t===1)return[Te(Ct(n[0],ut(-1)),n[1])];const e=[];let i=n.map(s=>ut(s.real,s.imag)),r=t*5;for(;i.length>2&&r-- >0;){const s=ut(.3+Math.random()*.7,.3+Math.random()*.7),{root:o}=$i(i,s,100),a=$i(n,o,20);e.push(a.root);const l=go(i,o);if(l.length>=i.length){console.warn("polyDeflate did not reduce degree, breaking");break}i=l}if(i.length===2)e.push(Te(Ct(i[0],ut(-1)),i[1]));else if(i.length===3){const s=i[2],o=i[1],a=i[0],l=dn(Ct(o,o),Ct(Ct(ut(4),s),a)),u=Rs(l),h=Ct(ut(2),s),c=Te(dn(Ct(ut(-1),o),u),h),m=Te(mn(Ct(ut(-1),o),u),h);e.push(c,m)}return e}function xo(n,t,e){const i=[ut(n),ut(-1),ut(n*t),ut(0),ut(n*e)],r=yo(i);if(r.length===0)return console.warn("laguerreSmallestPositiveRoot: no roots found"),n;let s=1/0,o=!1;for(const l of r)Math.abs(l.imag)<1e-10&&l.real>0&&l.real<s&&(s=l.real,o=!0);return o?$i(i,ut(s,0),20).root.real:(console.warn("laguerreSmallestPositiveRoot: no positive real root found"),n)}function bo(n,t,e){if(Math.abs(t)<1e-10&&Math.abs(e)<1e-10)return n;if(n<1e-10)return 0;if(Math.abs(e)<1e-10){const i=-1/(t*n),r=1/t,s=i*i-4*r;if(s<0)return n;const o=Math.sqrt(s),a=-.5*(i+Math.sign(i)*o),l=a,u=r/a;return l>0&&u>0?Math.min(l,u):l>0?l:u>0?u:n}try{return xo(n,t,e)}catch(i){return console.error("laguerreSmallestPositiveRoot failed:",i),n}}function _o(n,t,e,i,r){const s=n.x-t.x,o=n.y-t.y,a=Math.hypot(s,o)/Math.max(1e-12,e);if(a<1e-12)return{x:n.x,y:n.y};const l=a*a,u=1+(i+r*l)*l;return{x:s/u+t.x,y:o/u+t.y}}function wo(n,t,e,i,r){const s=n.x-t.x,o=n.y-t.y,a=Math.hypot(s,o)/Math.max(1e-12,e);if(a<1e-12)return{x:t.x,y:t.y};const u=bo(a,i,r)/a;return{x:t.x+s*u,y:t.y+o*u}}function Nt(n,t){const e=_o(n,{x:t.principalX,y:t.principalY},t.radiusNorm,t.k1,t.k2);return{x:e.x+(t.correctedOffsetX??0),y:e.y+(t.correctedOffsetY??0)}}function we(n,t){const e={x:n.x-(t.correctedOffsetX??0),y:n.y-(t.correctedOffsetY??0)};return wo(e,{x:t.principalX,y:t.principalY},t.radiusNorm,t.k1,t.k2)}const Mo=`
attribute vec2 a_position;

void main() {
  gl_Position = vec4(a_position, 0.0, 1.0);
}
`,So=`
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
`,Po=`
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
`;class vo{constructor(){yt(this,"canvas",null);yt(this,"gl",null);yt(this,"blurProgram",null);yt(this,"sobelProgram",null);yt(this,"positionBuffer",null);yt(this,"blurUniforms",null);yt(this,"sobelUniforms",null);yt(this,"resources",null);yt(this,"initialized",!1);yt(this,"unavailable",!1);yt(this,"maxTextureSize",0)}compute(t,e,i){if(!this.initialized&&!this.init())return null;const r=this.gl,s=this.blurProgram,o=this.sobelProgram,a=this.blurUniforms,l=this.sobelUniforms;if(!r||!s||!o||!a||!l||!this.positionBuffer||!this.canvas||e<=2||i<=2||e>this.maxTextureSize||i>this.maxTextureSize)return null;const u=this.ensureResources(e,i);if(!u)return null;this.canvas.width=e,this.canvas.height=i,r.viewport(0,0,e,i),r.disable(r.BLEND),r.pixelStorei(r.UNPACK_ALIGNMENT,1),r.pixelStorei(r.PACK_ALIGNMENT,1),r.bindBuffer(r.ARRAY_BUFFER,this.positionBuffer),r.activeTexture(r.TEXTURE0),r.bindTexture(r.TEXTURE_2D,u.sourceTexture),r.texImage2D(r.TEXTURE_2D,0,r.LUMINANCE,e,i,0,r.LUMINANCE,r.UNSIGNED_BYTE,t),r.useProgram(s),r.enableVertexAttribArray(0),r.vertexAttribPointer(0,2,r.FLOAT,!1,0,0),r.uniform2f(a.size,e,i),r.uniform1i(a.source,0),r.bindFramebuffer(r.FRAMEBUFFER,u.blurFramebuffer),r.bindTexture(r.TEXTURE_2D,u.sourceTexture),r.drawArrays(r.TRIANGLES,0,6);const h=new Uint8Array(e*i*4);r.readPixels(0,0,e,i,r.RGBA,r.UNSIGNED_BYTE,h),r.useProgram(o),r.uniform2f(l.size,e,i),r.uniform1i(l.blurred,0),r.bindFramebuffer(r.FRAMEBUFFER,u.sobelFramebuffer),r.bindTexture(r.TEXTURE_2D,u.blurTexture),r.drawArrays(r.TRIANGLES,0,6);const c=new Uint8Array(e*i*4);r.readPixels(0,0,e,i,r.RGBA,r.UNSIGNED_BYTE,c),r.disableVertexAttribArray(0),r.bindFramebuffer(r.FRAMEBUFFER,null),r.bindBuffer(r.ARRAY_BUFFER,null),r.bindTexture(r.TEXTURE_2D,null);const m=new Uint8Array(e*i);for(let y=0,x=0;y<m.length;y++,x+=4)m[y]=h[x];const d=new Float32Array(e*i),p=new Float32Array(e*i),f=new Float32Array(e*i);for(let y=0,x=0;y<d.length;y++,x+=4){const g=(c[x]|c[x+1]<<8)-32768,b=(c[x+2]|c[x+3]<<8)-32768;d[y]=g,p[y]=b,f[y]=Math.sqrt(g*g+b*b)}return{blurredGray:m,gx:d,gy:p,magnitude:f}}init(){if(this.initialized&&this.gl&&this.blurProgram&&this.sobelProgram)return!0;if(this.unavailable)return!1;const t=this.createCanvas();if(!t)return this.unavailable=!0,!1;const e=t.getContext("webgl",{alpha:!1,antialias:!1,depth:!1,stencil:!1,premultipliedAlpha:!1,preserveDrawingBuffer:!1});if(!e)return this.unavailable=!0,!1;const i=this.compileShader(e,e.VERTEX_SHADER,Mo),r=this.compileShader(e,e.FRAGMENT_SHADER,So),s=this.compileShader(e,e.FRAGMENT_SHADER,Po);if(!i||!r||!s)return i&&e.deleteShader(i),r&&e.deleteShader(r),s&&e.deleteShader(s),this.unavailable=!0,!1;const o=this.createProgram(e,i,r),a=this.createProgram(e,i,s);if(e.deleteShader(i),e.deleteShader(r),e.deleteShader(s),!o||!a)return o&&e.deleteProgram(o),a&&e.deleteProgram(a),this.unavailable=!0,!1;const l=e.createBuffer();return l?(e.bindBuffer(e.ARRAY_BUFFER,l),e.bufferData(e.ARRAY_BUFFER,new Float32Array([-1,-1,1,-1,-1,1,-1,1,1,-1,1,1]),e.STATIC_DRAW),e.bindBuffer(e.ARRAY_BUFFER,null),this.canvas=t,this.gl=e,this.blurProgram=o,this.sobelProgram=a,this.positionBuffer=l,this.blurUniforms={source:e.getUniformLocation(o,"u_source"),size:e.getUniformLocation(o,"u_size")},this.sobelUniforms={blurred:e.getUniformLocation(a,"u_blurred"),size:e.getUniformLocation(a,"u_size")},this.maxTextureSize=Number(e.getParameter(e.MAX_TEXTURE_SIZE)||0),this.initialized=!0,!0):(e.deleteProgram(o),e.deleteProgram(a),this.unavailable=!0,!1)}createCanvas(){return typeof OffscreenCanvas<"u"?new OffscreenCanvas(1,1):typeof document<"u"?document.createElement("canvas"):null}ensureResources(t,e){const i=this.gl;if(!i)return null;if(this.resources&&this.resources.width===t&&this.resources.height===e)return this.resources;this.disposeResources();const r=this.createTexture(i.LUMINANCE,t,e,i.LUMINANCE,i.UNSIGNED_BYTE,null),s=this.createTexture(i.RGBA,t,e,i.RGBA,i.UNSIGNED_BYTE,null),o=this.createTexture(i.RGBA,t,e,i.RGBA,i.UNSIGNED_BYTE,null),a=this.createFramebuffer(s),l=this.createFramebuffer(o);return!r||!s||!o||!a||!l?(r&&i.deleteTexture(r),s&&i.deleteTexture(s),o&&i.deleteTexture(o),a&&i.deleteFramebuffer(a),l&&i.deleteFramebuffer(l),null):(this.resources={width:t,height:e,sourceTexture:r,blurTexture:s,sobelTexture:o,blurFramebuffer:a,sobelFramebuffer:l},this.resources)}createTexture(t,e,i,r,s,o){const a=this.gl;if(!a)return null;const l=a.createTexture();return l?(a.bindTexture(a.TEXTURE_2D,l),a.texParameteri(a.TEXTURE_2D,a.TEXTURE_MIN_FILTER,a.NEAREST),a.texParameteri(a.TEXTURE_2D,a.TEXTURE_MAG_FILTER,a.NEAREST),a.texParameteri(a.TEXTURE_2D,a.TEXTURE_WRAP_S,a.CLAMP_TO_EDGE),a.texParameteri(a.TEXTURE_2D,a.TEXTURE_WRAP_T,a.CLAMP_TO_EDGE),a.texImage2D(a.TEXTURE_2D,0,t,e,i,0,r,s,o),a.bindTexture(a.TEXTURE_2D,null),l):null}createFramebuffer(t){const e=this.gl;if(!e||!t)return null;const i=e.createFramebuffer();if(!i)return null;e.bindFramebuffer(e.FRAMEBUFFER,i),e.framebufferTexture2D(e.FRAMEBUFFER,e.COLOR_ATTACHMENT0,e.TEXTURE_2D,t,0);const r=e.checkFramebufferStatus(e.FRAMEBUFFER);return e.bindFramebuffer(e.FRAMEBUFFER,null),r!==e.FRAMEBUFFER_COMPLETE?(e.deleteFramebuffer(i),null):i}compileShader(t,e,i){const r=t.createShader(e);return r?(t.shaderSource(r,i),t.compileShader(r),t.getShaderParameter(r,t.COMPILE_STATUS)?r:(console.error("[SFR Auto Detect WebGL] shader compile failed",t.getShaderInfoLog(r)),t.deleteShader(r),null)):null}createProgram(t,e,i){const r=t.createProgram();return r?(t.attachShader(r,e),t.attachShader(r,i),t.bindAttribLocation(r,0,"a_position"),t.linkProgram(r),t.getProgramParameter(r,t.LINK_STATUS)?r:(console.error("[SFR Auto Detect WebGL] program link failed",t.getProgramInfoLog(r)),t.deleteProgram(r),null)):null}disposeResources(){const t=this.gl,e=this.resources;if(!t||!e){this.resources=null;return}t.deleteTexture(e.sourceTexture),t.deleteTexture(e.blurTexture),t.deleteTexture(e.sobelTexture),t.deleteFramebuffer(e.blurFramebuffer),t.deleteFramebuffer(e.sobelFramebuffer),this.resources=null}}const Co=new vo,Fo=12,ko=12,Ao=501,To=2,Io=[.3,.4,.55,.75,1,1.35],No=320,Ro=.75,Lo=.002,Ls=.75,Es=2,ai=2,ts=32,St=1e-12;function Ft(n,t,e){return Math.max(t,Math.min(e,n))}function fn(n){if(n.length===0)return NaN;const t=[...n].sort((i,r)=>i-r),e=t.length>>1;return t.length&1?t[e]:.5*(t[e-1]+t[e])}function cr(n){if(n.length===0)return 0;const t=fn(n),e=n.map(i=>Math.abs(i-t));return 1.4826*fn(e)}function Eo(n,t){const e=n.map(a=>(a.x%1+1)%1).sort((a,l)=>a-l),i=[];for(const a of e)(i.length===0||a-i[i.length-1]>1e-5)&&i.push(a);let r=1;if(i.length>1){r=1-i[i.length-1]+i[0];for(let a=1;a<i.length;a++)r=Math.max(r,i[a]-i[a-1])}const s=Ft(1-r,0,1);let o=t;return i.length<8||s<.75?o=Math.min(o,.5):s<.9&&(o=Math.min(o,1)),n.length<96&&(o=Math.min(o,.75)),{uniquePhaseCount:i.length,maximumPhaseGapPx:r,phaseCoverage:s,reliableFrequency:o}}function Uo(n,t){let e=0,i=n.length;for(;e<i;){const r=e+i>>1;n[r].x<t?e=r+1:i=r}return e}function Do(n,t){let e=0,i=n.length;for(;e<i;){const r=e+i>>1;n[r].x<=t?e=r+1:i=r}return e}function ur(n,t){const e=t.length,i=n.map((r,s)=>[...r,t[s]]);for(let r=0;r<e;r++){let s=r;for(let a=r+1;a<e;a++)Math.abs(i[a][r])>Math.abs(i[s][r])&&(s=a);if(Math.abs(i[s][r])<=1e-14)return null;s!==r&&([i[s],i[r]]=[i[r],i[s]]);const o=i[r][r];for(let a=r;a<=e;a++)i[r][a]/=o;for(let a=0;a<e;a++){if(a===r)continue;const l=i[a][r];if(l!==0)for(let u=r;u<=e;u++)i[a][u]-=l*i[r][u]}}return i.map(r=>r[e])}function Bo(n,t,e,i,r,s=-1){const o=Uo(n,t-e),a=Do(n,t+e),l=i+1,u=Array.from({length:l},()=>new Array(l).fill(0)),h=new Array(l).fill(0);let c=0,m=0;for(let f=o;f<a;f++){const y=n[f];if(s>=0&&y.fold===s)continue;const x=(y.x-t)/e,g=Math.abs(x);if(g>=1)continue;const _=Math.pow(1-g*g*g,3)*(r?r[f]:1);if(!(_>0)||!Number.isFinite(_))continue;const M=new Array(l).fill(1);for(let w=1;w<l;w++)M[w]=M[w-1]*x;for(let w=0;w<l;w++){h[w]+=_*M[w]*y.y;for(let P=0;P<l;P++)u[w][P]+=_*M[w]*M[P]}c++,m+=_}if(c<Math.max(8,l*2)||m<=St)return null;const d=Math.max(St,u[0][0]*1e-11);for(let f=0;f<l;f++)u[f][f]+=d*Math.pow(4,f);const p=ur(u,h);return!p||p.some(f=>!Number.isFinite(f))?null:{value:p[0],derivative:i>=1?p[1]/e:0,effectiveCount:c,weightSum:m}}function Oo(n,t){const e=t*.65,i=n.filter(o=>Math.abs(o.x)>=e);if(i.length<12)return 0;let r=new Array(i.length).fill(1),s=null;for(let o=0;o<3;o++){const a=Array.from({length:3},()=>new Array(3).fill(0)),l=new Array(3).fill(0);for(let c=0;c<i.length;c++){const m=i[c],d=[1,m.x>=0?1:0,m.x/t],p=r[c];for(let f=0;f<3;f++){l[f]+=p*d[f]*m.y;for(let y=0;y<3;y++)a[f][y]+=p*d[f]*d[y]}}if(s=ur(a,l),!s)return 0;const u=i.map(c=>{const m=s[0]+s[1]*(c.x>=0?1:0)+s[2]*c.x/t;return c.y-m}),h=Math.max(St,cr(u));r=u.map(c=>{const m=Math.abs(c)/(1.345*h);return m<=1?1:1/m})}return s?s[2]/t:0}function zo(n,t,e,i){const r=[];for(let g=0;g<n.length;g++){const b=Number(n[g]),_=Number(t[g]);Number.isFinite(b)&&Number.isFinite(_)&&r.push({x:b,y:_,fold:g%5})}if(r.sort((g,b)=>g.x-b.x),r.length<ts)return null;const s=Math.min(Math.abs(r[0].x),Math.abs(r[r.length-1].x)),o=Ft(Number.isFinite(e)?e:Math.min(ko,s*.96),3,Math.max(3,s*.98));let a=r.filter(g=>Math.abs(g.x)<=o);if(a.length<ts)return null;const l=Ft(Math.floor(a.length*.08),8,64),u=fn(a.slice(0,l).map(g=>g.y)),h=fn(a.slice(-l).map(g=>g.y));if(!Number.isFinite(u)||!Number.isFinite(h)||Math.abs(h-u)<=St)return null;u>h&&(a=a.map(g=>({...g,x:-g.x})).sort((g,b)=>g.x-b.x));const c=i?Oo(a,o):0,m=a.map(g=>({...g,y:g.y-c*g.x})),d=fn(m.slice(0,l).map(g=>g.y)),f=fn(m.slice(-l).map(g=>g.y))-d;if(!Number.isFinite(f)||f<=St)return null;const y=m.map(g=>({...g,y:(g.y-d)/f})),x=[...y.slice(0,l).map(g=>g.y),...y.slice(-l).map(g=>g.y-1)];return{samples:y,halfWidth:o,normalization:{dark:d,contrast:f,removedSlopePerPixel:c,noiseSigmaNormalized:cr(x)}}}function hr(n,t,e){if(n<=t[0])return e[0];if(n>=t[t.length-1])return e[e.length-1];let i=0,r=t.length-1;for(;r-i>1;){const o=i+r>>1;t[o]<=n?i=o:r=o}const s=(n-t[i])/(t[r]-t[i]);return e[i]+s*(e[r]-e[i])}function Nn(n,t,e,i,r,s=null,o=t){const a=1/e,l=Ft(o,Math.min(Es,t),t),u=Math.max(33,Math.round(2*l/a)+1),h=new Array(u);for(let x=0;x<u;x++)h[x]=-l+x*(2*l/(u-1));const c=s?n.map((x,g)=>({...x,y:s(x,g)})):n;let m=new Array(c.length).fill(1),d=[],p=[],f=[];const y=()=>{d=new Array(u),p=new Array(u),f=new Array(u);for(let x=0;x<u;x++){let g=i,b=null;for(let _=0;_<4&&!b;_++)b=Bo(c,h[x],g,3,m),g*=1.35;if(!b)return!1;d[x]=b.value,p[x]=b.derivative,f[x]=b.effectiveCount}return!0};if(!y())return null;if(!s)for(let x=0;x<r;x++){const g=c.map(_=>_.y-hr(_.x,h,d)),b=Math.max(St,cr(g));if(m=g.map(_=>{const M=_/(4.685*b);if(Math.abs(M)>=1)return .001;const w=1-M*M;return Math.max(.001,w*w)}),!y())return null}return{x:h,esf:d,lsf:p,supportCounts:f,step:h[1]-h[0],profileHalfWidth:l}}function Vo(n,t,e){const i=n.x.reduce((p,f,y)=>Math.abs(f)<Math.abs(n.x[p])?y:p,0),r=Ft(4*t,.035,.12),s=Math.max(1.25,1.5*e),o=(p,f)=>{let y=Math.abs(n.x[i]),x=0,g=!1;for(let b=i;b>=0&&b<n.x.length;b+=p)if(Math.abs(n.esf[b]-f)>r)y=Math.abs(n.x[b]),x=0;else if(x+=n.step,x>=s){g=!0;break}return{extentPx:y,foundStablePlateau:g}},a=o(-1,0),l=o(1,1),u=Math.max(a.extentPx,l.extentPx),h=Math.max(3,4*e),c=n.profileHalfWidth*.96,m=Ft(u+h,Math.min(ai,c),c),d=Math.min(u+.5*h,m*.82);return{transitionHalfWidthPx:u,windowHalfWidthPx:m,flatFraction:Ft(d/Math.max(St,m),.25,.82),plateauTolerance:r,stablePlateauRequiredPx:s,boundaryLimited:!a.foundStablePlateau||!l.foundStablePlateau}}function Go(n){let t=1;for(;t<n;)t*=2;return t}function Xo(n,t=null){const e=n.length;if(e<1||e&e-1)throw new Error("FFT length must be a power of two.");const i=Float64Array.from(n),r=t?Float64Array.from(t):new Float64Array(e);for(let s=1,o=0;s<e;s++){let a=e>>1;for(;o&a;)o^=a,a>>=1;o^=a,s<o&&([i[s],i[o]]=[i[o],i[s]],[r[s],r[o]]=[r[o],r[s]])}for(let s=2;s<=e;s*=2){const o=-2*Math.PI/s,a=Math.cos(o),l=Math.sin(o);for(let u=0;u<e;u+=s){let h=1,c=0;const m=s>>1;for(let d=0;d<m;d++){const p=u+d,f=p+m,y=i[f]*h-r[f]*c,x=i[f]*c+r[f]*h,g=i[p],b=r[p];i[p]=g+y,r[p]=b+x,i[f]=g-y,r[f]=b-x;const _=h*a-c*l;c=h*l+c*a,h=_}}}return{real:i,imag:r}}function Yo(n,t,e,i){const r=Math.abs(n-t)/Math.max(St,e);if(r>=1)return 0;if(r<=i)return 1;const s=(r-i)/(1-i);return .5*(1+Math.cos(Math.PI*s))}function oi(n,t,e=null,i=null,r=null,s=.72){var C;const o=Number.isFinite(i)?i:0,a=Math.min(o-n.x[0],n.x[n.x.length-1]-o),l=Number.isFinite(r)?Ft(r,Math.min(ai,a),a):Math.max(ai,a*.98),u=n.x.map((v,k)=>({position:v,index:k})).filter(v=>Math.abs(v.position-o)<=l),h=u.reduce((v,k)=>Math.abs(n.lsf[k.index])>Math.abs(n.lsf[v])?k.index:v,((C=u[0])==null?void 0:C.index)??0),c=n.lsf.reduce((v,k,A)=>Math.abs(k)>Math.abs(n.lsf[v])?A:v,0),m=n.x[h],d=n.x[c],p=Ft(s,0,.9),f=n.lsf.map((v,k)=>v*Yo(n.x[k],o,l,p)),y=e||Go(Math.max(2048,f.length*8)),x=new Float64Array(y);x.set(f.slice(0,y));const g=Xo(x),b=Math.max(St,Math.hypot(g.real[0],g.imag[0])),_=1/(y*n.step),M=Math.min(Math.floor(t/_)+1,y>>1),w=new Array(M),P=new Array(M);for(let v=0;v<M;v++)w[v]=v*_,P[v]=Math.hypot(g.real[v],g.imag[v])/b;return{frequencies:w,mtf:P,fftSize:y,center:o,peakPosition:m,globalPeakPosition:d,windowHalfWidth:l,flatFraction:p}}function Wo(n,t,e,i,r,s,o,a,l,u,h,c,m=new Map){const d=i.map(Number).filter(_=>Number.isFinite(_)&&_>=r-St&&_<=Math.max(Ls,r)+St).sort((_,M)=>_-M);d.some(_=>Math.abs(_-r)<=St)||d.unshift(r);const p=d.filter((_,M)=>M===0||Math.abs(_-d[M-1])>St),f=Math.min(l,.3),y=[];for(let _=.05;_<=f+St;_+=.05)y.push(_);const x=[];for(const _ of p){const M=_.toFixed(8),w=m.get(M)??Nn(n,t,e,_,s,null,o);w&&m.set(M,w);const P=w?Nn(n,t,e,_,0,A=>A.x>=0?1:0,o):null;if(!w||!P){x.push({bandwidth:_,profile:w,idealProfile:P,curve:null,peakOffsetPx:1/0});continue}const C=oi(w,f,null,0,a.windowHalfWidthPx,a.flatFraction),v=oi(P,f,C.fftSize,0,a.windowHalfWidthPx,a.flatFraction),k=C.mtf.map((A,S)=>{const F=Math.max(u,v.mtf[S]||u);return A*Math.min(h,1/F)});x.push({bandwidth:_,profile:w,idealProfile:P,curve:y.map(A=>hr(A,C.frequencies,k)),peakOffsetPx:Math.abs(C.peakPosition)})}const g=Ft(Math.max(.015,2.25*c),.015,.05);for(let _=0;_<x.length;_++){const M=x[_];let w=M.curve&&M.peakOffsetPx<=2?0:1/0;if(Number.isFinite(w))for(let P=_+1;P<x.length;P++){const C=x[P];if(!C.curve||C.peakOffsetPx>2){w=1/0;break}for(let v=0;v<M.curve.length;v++)w=Math.max(w,Math.abs(M.curve[v]-C.curve[v]))}if(M.maximumCoarserDifference=w,w<=g)return{bandwidth:M.bandwidth,profile:M.profile,idealProfile:M.idealProfile,tolerance:g,stabilityMaximumFrequency:f,evaluations:x}}const b=x[x.length-1];return{bandwidth:b.bandwidth,profile:b.profile,idealProfile:b.idealProfile,tolerance:g,stabilityMaximumFrequency:f,evaluations:x}}function Ho(n,t,e){for(let i=1;i<t.length;i++)if(t[i-1]>e&&t[i]<=e){const r=t[i]-t[i-1],s=Math.abs(r)<=St?1:(e-t[i-1])/r;return{frequency:n[i-1]+s*(n[i]-n[i-1]),index:i}}return null}function jo(n,t,e,i){if(!e)return null;const r=n[1]-n[0],s=Math.max(8,Math.round(.035/Math.max(St,r))),o=Math.max(0,e.index-s),a=Math.min(t.length-1,e.index+s);if(a-o+1<9)return null;const l=e.frequency,u=Math.max(r,n[a]-n[o]),h=Array.from({length:3},()=>new Array(3).fill(0)),c=new Array(3).fill(0);for(let g=o;g<=a;g++){if(!(t[g]>0)||!Number.isFinite(t[g]))continue;const b=(n[g]-l)/u,_=Math.exp(-4*b*b),M=[1,b,b*b],w=Math.log(t[g]);for(let P=0;P<3;P++){c[P]+=_*M[P]*w;for(let C=0;C<3;C++)h[P][C]+=_*M[P]*M[C]}}const m=ur(h,c);if(!m||m[1]>=0)return null;const d=m[0]-Math.log(i),p=m[1],f=m[2],y=[];if(Math.abs(f)<=1e-10)Math.abs(p)>St&&y.push(-d/p);else{const g=p*p-4*f*d;if(g>=0){const b=Math.sqrt(g);y.push((-p-b)/(2*f),(-p+b)/(2*f))}}return y.map(g=>l+g*u).filter(g=>Number.isFinite(g)&&g>=n[o]&&g<=n[a]).sort((g,b)=>Math.abs(g-l)-Math.abs(b-l))[0]??null}function Xi(n,t,e){return n.map(i=>hr(i,t,e))}function qo(n,t,e){for(let i=1;i<t.length;i++)if(t[i-1]>=e&&t[i]<e){const r=t[i]-t[i-1],s=Math.abs(r)<=St?0:(e-t[i-1])/r;return n[i-1]+s*(n[i]-n[i-1])}return n[n.length-1]}function Us(n,t,e={}){if(!n||!t||n.length!==t.length)throw new TypeError("distancesPx and values must have the same length.");const i=Ft(Math.round(Number(e.oversampling)||Fo),4,16),r=Ft(Number(e.maxFrequencyCyclesPerPixel)||To,.5,i*.45),s=Ft(Math.round(Number(e.outputPoints)||Ao),65,4001),o=e.detrend!==!1,a=Number(e.halfWidthPx),l=zo(n,t,a,o);if(!l)return null;const u=l.samples.length<No,h=Number(e.robustIterations),c=Number.isFinite(h)?Ft(Math.round(h),0,4):l.normalization.noiseSigmaNormalized>=Lo?u?4:2:0,m=Ft(Number(e.correctionFloor)||.18,.05,.8),d=Ft(Number(e.correctionMaxGain)||3,1,8),p=(e.bandwidthCandidatesPx||Io).map(Number).filter(L=>Number.isFinite(L)&&L>=.2&&L<=2.5).sort((L,G)=>L-G),f=Number(e.bandwidthPx),y=Number.isFinite(f)?Ft(f,.2,2.5):u?Ro:p[0]??.55,x=Ls,g=Math.max(x*.55,2/i),b=l.halfWidth-g;if(b<Es)return null;const _=Nn(l.samples,l.halfWidth,i,x,c,null,b);if(!_)return null;const M=Vo(_,l.normalization.noiseSigmaNormalized,x);if(M.windowHalfWidthPx<ai)return null;const w=new Map([[x.toFixed(8),_]]);let P=y;const C=Number.isFinite(f)?"requested":u?"small-window-fixed":"multiscale-derivative-stability";let v=null,k=w.get(P.toFixed(8))??null,A=null;if(!Number.isFinite(f)&&!u&&(v=Wo(l.samples,l.halfWidth,i,p,y,c,b,M,r,m,d,l.normalization.noiseSigmaNormalized,w),P=v.bandwidth,k=v.profile,A=v.idealProfile),k||(k=Nn(l.samples,l.halfWidth,i,P,c,null,b)),!k||(A=A??Nn(l.samples,l.halfWidth,i,P,0,L=>L.x>=0?1:0,b),!A))return null;const S=oi(k,r,null,0,M.windowHalfWidthPx,M.flatFraction);if(Math.abs(S.peakPosition)>2)return null;const F=oi(A,r,S.fftSize,0,M.windowHalfWidthPx,M.flatFraction),T=S.mtf.map((L,G)=>{const X=Math.max(m,F.mtf[G]||m);return L*Math.min(d,1/X)}),R=Eo(l.samples,r),E=qo(F.frequencies,F.mtf,1/d),O=Math.min(R.reliableFrequency,E),U=Ft(Number(e.mtfLevel)||.5,.05,.95),I=Ho(S.frequencies,T,U),D=jo(S.frequencies,T,I,U),B=new Array(s);for(let L=0;L<s;L++)B[L]=L*r/(s-1);const V=Xi(B,S.frequencies,S.mtf),Q=Xi(B,S.frequencies,T),z=Xi(B,F.frequencies,F.mtf),Y=Number(e.pixelPitchUm),W=Number.isFinite(Y)&&Y>0,K=W?1e3/Y:1,it=B.map(L=>L*K),J=D??(I==null?void 0:I.frequency)??null,ot=Math.min(...z.filter(Number.isFinite)),Z=[];return Math.abs(S.peakPosition)>1&&Z.push("LSF peak is more than one pixel from the supplied edge origin."),M.boundaryLimited&&Z.push("A stable plateau was not established on both sides before the profile boundary; low-frequency confidence is reduced."),ot<m&&Z.push("High-frequency estimator correction reached its configured floor."),R.phaseCoverage<.75?Z.push("Projected samples do not cover a complete subpixel phase cycle; frequencies above 0.5 cycles/pixel are not reliable."):R.phaseCoverage<.9&&Z.push("Projected samples have incomplete subpixel phase coverage; high-frequency results have reduced confidence."),E<r-St&&Z.push("Estimator correction reached its configured gain limit before the requested maximum frequency."),J!==null&&J>O+St&&Z.push("The requested MTF crossing is above the reliable frequency limit for this sample geometry."),I||Z.push(`MTF did not cross ${U.toFixed(3)} inside the requested frequency range.`),{esf:k.esf,lsf:k.lsf,profilePositionsPx:k.x,frequencies:it,mtf:Q,mtfRaw:V,estimatorResponse:z,mtf50:J===null?null:J*K,mtf50Linear:I?I.frequency*K:null,mtf50Method:D===null?"linear":"local-log-quadratic",frequencyUnit:W?"cycles/mm":"cycles/pixel",diagnostics:{sampleCount:l.samples.length,oversampling:i,halfWidthPx:l.halfWidth,profileHalfWidthPx:k.profileHalfWidth,boundaryGuardPx:g,bandwidthPx:P,preliminaryBandwidthPx:y,bandwidthSelectionMethod:C,bandwidthStabilityTolerance:(v==null?void 0:v.tolerance)??null,bandwidthStabilityMaximumFrequencyCyclesPerPixel:(v==null?void 0:v.stabilityMaximumFrequency)??null,bandwidthStabilityCandidates:(v==null?void 0:v.evaluations.map(L=>({bandwidthPx:L.bandwidth,peakOffsetPx:L.peakOffsetPx,maximumCoarserDifference:L.maximumCoarserDifference??null})))??[],pilotBandwidthPx:x,smallWindowMode:u,robustIterations:c,profileStepPx:k.step,fftSize:S.fftSize,lsfPeakPositionPx:S.peakPosition,lsfGlobalPeakPositionPx:S.globalPeakPosition,windowHalfWidthPx:S.windowHalfWidth,windowFlatFraction:S.flatFraction,transitionHalfWidthPx:M.transitionHalfWidthPx,plateauTolerance:M.plateauTolerance,stablePlateauRequiredPx:M.stablePlateauRequiredPx,plateauBoundaryLimited:M.boundaryLimited,correctionFloor:m,correctionMaxGain:d,samplingUniquePhaseCount:R.uniquePhaseCount,samplingMaximumPhaseGapPx:R.maximumPhaseGapPx,samplingPhaseCoverage:R.phaseCoverage,estimatorReliableFrequencyCyclesPerPixel:E,maxReliableFrequencyCyclesPerPixel:O,normalization:l.normalization,minLocalSupport:Math.min(...k.supportCounts),warnings:Z}}}const un=(n,t)=>{const e=$e(n);if(!e)throw new Error(`Bayer CFA pattern is unresolved for ${t}.`);return e},Ji={gradientPercentiles:[.82,.88,.92,.95,.98,.995],downsampleMaxSide:1600,minComponentAreaRatio:15e-6,maxComponentAreaRatio:.35,minComponentAreaPx:20,minEdgePoints:24,extentQuantileLow:.02,extentQuantileHigh:.98,cornerTrimRatio:.18,minSpanPx:8,maxAspectRatio:2,bandScale:.16,bandMinPx:1.75,bandMaxPx:14,minPointContrast:6,minSidePoints:3,minCoverageRatio:.15,minCenterCoverageRatio:.2,filterBlockPurity:!0,innerPurityStdScale:1.5,outerMeanSpreadLimit:51,minAxisDot:.6,residualLimitFloor:.01,residualLimitScale:.25,minQuadArea:48,minSideLength:10,minOuterContrast:5,sampleHalfWidthRatio:.25};function Qo(n,t,e,i,r,s){const o=n.width,a=n.height,l=un(n.bayerPattern,"corrected RAW SFR sampling"),u=[],h=[],c=s!=null&&s.correctedRect?Tt*2:Tt,m=Math.max(1,Math.min(r,c)),d=e.p2.x-e.p1.x,p=e.p2.y-e.p1.y,f=Math.hypot(d,p);if(!Number.isFinite(f)||f<=1e-6)return null;const y=d/f,x=p/f,g=-x,b=y,_={x:(e.p1.x+e.p2.x)*.5,y:(e.p1.y+e.p2.y)*.5},M=s!=null&&s.correctedRect?xt(s.correctedRect,o,a):xt(Et(de(e,c*4+2)??[e.p1,e.p2],2),o,a);if(!M)return null;const w=(s==null?void 0:s.correctedScanlinesOverride)??(s!=null&&s.distortedRect?_r(xt(s.distortedRect,o,a)??s.distortedRect,t,o,a):$s(M,e,Math.max(1,i),m*4+.5,o,a));if(!w||w.size===0)return null;const P=Js(w,t,o,a);if(P.size===0)return null;const C=!wr(t);for(const[k,A]of P){if(k<0||k>=a)continue;const S=k*o;for(let F=A.start;F<=A.end;F++){if(F<0||F>=o||!_t(F,k,l,s==null?void 0:s.greenPhase))continue;const T={x:F,y:k},R=Nt(T,t);if(!Number.isFinite(R.x)||!Number.isFinite(R.y)||Math.round(R.x)<0||Math.round(R.x)>=o||Math.round(R.y)<0||Math.round(R.y)>=a)continue;const E=R.x-_.x,O=R.y-_.y,U=E*y+O*x;let I=E*g+O*b;if(C){const D=Zs(U,y,x,_,T,t);if(!D)continue;const B=.5*(D.a+D.b),V=Yt(D.a,y,x,_,t),Q=Yt(B,y,x,_,t),z=Yt(D.b,y,x,_,t),Y=Sr({x:D.a,y:Math.hypot(V.x-T.x,V.y-T.y)},{x:B,y:Math.hypot(Q.x-T.x,Q.y-T.y)},{x:D.b,y:Math.hypot(z.x-T.x,z.y-T.y)});if(!Number.isFinite(Y))continue;const W=Mr(Y,y,x,_,t),K=Math.hypot(W.x,W.y);if(!Number.isFinite(K)||K<=1e-9)continue;const it=W.x/K,ot=-(W.y/K),Z=it,L=Yt(Y,y,x,_,t);I=(T.x-L.x)*ot+(T.y-L.y)*Z}!Number.isFinite(U)||Math.abs(U)>Math.max(1,i)||!Number.isFinite(I)||Math.abs(I)>m||(u.push(I),h.push(Math.max(0,n.data[S+F]-bn(s==null?void 0:s.blackLevel,F,k))))}}if(u.length<8)return null;const v=Math.abs(d)>=Math.abs(p)?1:2;return s!=null&&s.forceLegacyModel?wn(u,h,v,c):_n(u,h,v,c)}function dr(n,t){const e=n.length;let i=0,r=0,s=0,o=0;for(let l=0;l<e;l++)i+=n[l],r+=t[l],s+=n[l]*t[l],o+=n[l]*n[l];const a=e*o-i*i;return a===0?{slope:0,intercept:0}:{slope:(e*s-i*r)/a,intercept:(r*o-i*s)/a}}function Ko(n,t){const e=n.length,i=new Array(e).fill(0),r=2*t;for(let s=0;s<e;s++){const o=s>0?n[s-1]:n[0],a=s<e-1?n[s+1]:n[e-1];i[s]=(a-o)/r}return i}const An=-1e7,li=13,ue=512,It=8,fr=1/It,Tt=28,$o=[[0,0,0,0,0,-.085714285714286,.342857142857143,.485714285714286,.342857142857143,-.085714285714286,0,0,0,0,0],[0,0,0,0,-.095238095238095,.142857142857143,.285714285714286,.333333333333333,.285714285714286,.142857142857143,-.095238095238095,0,0,0,0],[0,0,0,-.090909090909091,.060606060606061,.168831168831169,.233766233766234,.255411255411255,.233766233766234,.168831168831169,.060606060606061,-.090909090909091,0,0,0],[0,0,-.083916083916084,.020979020979021,.102564102564103,.160839160839161,.195804195804196,.207459207459208,.195804195804196,.160839160839161,.102564102564103,.020979020979021,-.083916083916084,0,0],[0,-.076923076923077,0,.062937062937063,.111888111888112,.146853146853147,.167832167832168,.174825174825175,.167832167832168,.146853146853147,.111888111888112,.062937062937063,0,-.076923076923077,0],[-.070588235294118,-.011764705882353,.038009049773756,.078733031674208,.110407239819004,.133031674208145,.146606334841629,.151131221719457,.146606334841629,.133031674208145,.110407239819004,.078733031674208,.038009049773756,-.011764705882353,-.070588235294118]];function ci(n,t,e,i=1){const r=Math.max(1e-6,e*.5),s=Math.max(1e-6,t*i),o=Math.exp(-s*r),a=1-o;if(!Number.isFinite(a)||Math.abs(a)<=1e-9)return Math.abs(n)<=r?1:0;if(Math.abs(n)<r){const l=2-2*o*Math.cosh(s*n),u=2*Math.sinh(s*r)*a;return!Number.isFinite(u)||Math.abs(u)<=1e-9?0:l/u}return Math.exp(-s*Math.abs(n))/a}function Jo(n,t,e,i,r,s){const o=n.length;if(o===0)return[];if(s<1)return n;const a=Math.min(s,32),l=new Array(o).fill(0);l[0]=n[0];for(let d=1;d<o;d++)l[d]=l[d-1]+n[d];const u=(d,p)=>{const f=Math.max(0,d),y=Math.min(o-1,p);return y<f?n[Math.max(0,Math.min(o-1,d))]??0:(l[y]-(f>0?l[f-1]:0))/(y-f+1)},h=a*2,c=a,m=1;for(let d=Math.max(t+c,i-h);d<i;d++){const p=Math.max(m,Math.trunc((i-d)*c/Math.max(1,h)));n[d]=u(d-p,d+p)}for(let d=Math.min(r+h-1,e-c-1);d>r;d--){const p=Math.max(m,Math.trunc((d-r)*c/Math.max(1,h)));n[d]=u(d-p,d+p)}for(let d=c+1;d<i-h;d++)n[d]=u(d-c,d+c);for(let d=Math.min(r+h,e-c-1);d<o-c-1;d++)n[d]=u(d-c,d+c);return n}function Ds(n){return!Number.isFinite(n)||Math.abs(n)<=1e-9?1:Math.sin(n)/n}let Qn=null,Yi=null;function Zo(){if(Qn)return Qn;const n=.625,t=1/128,e=Math.max(16,Math.round(n*2/t)+1),i=[],r=[];for(let c=0;c<e;c++){const m=-n+c*t;i.push(m),r.push(Math.abs(m)<=n?ci(m,li,.125,1):0)}const s=4,o=1/1024,a=Math.round(s/o)+1,l=new Array(a).fill(0),u=new Array(a).fill(1);let h=0;for(let c=0;c<r.length;c++)h+=r[c];h=Math.max(1e-9,h);for(let c=0;c<a;c++){const m=c*o;l[c]=m;let d=0;for(let p=0;p<i.length;p++)d+=r[p]*Math.cos(2*Math.PI*m*i[p]);u[c]=Math.max(1e-6,Math.abs(d)/h)}return Qn={freqs:l,values:u},Qn}function tl(n,t){const e=Math.max(1e-6,Ds(Math.PI*n*t)),i=Zo(),r=Ie(Math.max(0,Math.min(i.freqs[i.freqs.length-1],n)),i.freqs,i.values);return Math.max(1e-6,e*r)}function el(){if(Yi)return Yi;const n=new Array(ue/16*4).fill(1),t=ue*16,e=new Float32Array(t);for(let s=0;s<t;s++){const o=(s-t/2)/(16*It);e[s]=Math.abs(o)<=.625?ci(o,li,fr,1):0}const i=new gn(t);i.transform(e);const r=Math.max(1e-9,Math.abs(i._real[0]));n[0]=1;for(let s=1;s<n.length;s++){const o=Ds(Math.PI*s/256),a=Math.max(1e-6,Math.hypot(i._real[s],i._imag[s])/r);n[s]=Math.max(1e-6,o*a)}return Yi=n,n}function nl(n,t,e=Tt){if(n.length===0||t.length!==n.length)return null;const i=ue,r=i/2,s=fr,o=2*It,a=Math.max(1,Math.round(.5*It)),l=5,u=Math.max(0,Math.round(r-e*It)),h=Math.min(i-1,Math.round(r+e*It));if(h-u<32)return null;let c=new Array(i).fill(An),m=0,d=0,p=0,f=0,y=u,x=h,g=u,b=h,_=0;for(;;){const L=new Array(i).fill(0),G=new Array(i).fill(0);c=new Array(i).fill(An),m=0,d=0,p=0,f=0;let X=-1,j=-1;for(let H=0;H<n.length;H++){const pt=Math.trunc(n[H]*It+r),mt=Math.max(g,pt-l),gt=Math.min(b-1,pt+l);for(let ht=mt;ht<=gt;ht++){const Zt=(ht-r)*s,Me=Math.max(0,1-Math.abs((n[H]-Zt)*1.75));Me>0&&(G[ht]+=t[H]*Me,L[ht]+=Me)}}const $=Math.max(r-i/8,g+2*It),et=Math.min(r+i/8,b-2*It);for(let H=Math.max(0,g-1);H<=Math.min(i-1,b+1);H++)L[H]>0&&(c[H]=G[H]/L[H],H<$&&(m+=c[H],p++),H>et&&(d+=c[H],f++),X<0&&(X=H),j=H);if(X<0||j<0||p<=0||f<=0)return null;for(let H=X-1;H>=0;H--)c[H]=c[X];for(let H=j+1;H<i;H++)c[H]=c[j];const nt=new Array(i).fill(0);let ct=r;const at=2*It;for(let H=at+1;H<i-1-at;H++){let pt=0,mt=0;for(let gt=-at;gt<=at;gt++)mt+=c[H+gt]*gt,pt+=gt*gt;nt[H]=pt>0?mt/pt:0,Math.abs(nt[H])>Math.abs(nt[ct]??0)&&H>g+at&&H<b-at-1&&(ct=H)}if(Math.abs(ct-r)>2*It&&Math.abs(ct-r)<12*It)return null;let rt=0;for(let H=Math.max(0,r-at);H<=Math.min(i-1,r+at);H++)Math.abs(nt[H])>Math.abs(rt)&&(rt=nt[H]);if(!Number.isFinite(rt)||Math.abs(rt)<=1e-9)return null;const Jt=Math.abs(rt*.001);y=g,x=b;let Wt=!1;for(let H=r-o;H>=g+a;H--)if(nt[H]*rt<0&&Math.abs(nt[H])>Jt){let pt=0,mt=0,gt=0;for(let ht=H;ht>=g;ht--)nt[ht]*rt<0&&(pt++,mt=Math.max(mt,Math.abs(nt[ht]))),gt++;if(pt>gt*.4&&mt/Math.abs(rt)>.25||pt>.9*gt&&gt>o){y=Math.min(H,r-o),Wt=!0;break}}for(let H=r+o;H<b-a;H++)if(nt[H]*rt<0&&Math.abs(nt[H])>Jt){let pt=0,mt=0,gt=0;for(let ht=H;ht<b;ht++)nt[ht]*rt<0&&(pt++,mt=Math.max(mt,Math.abs(nt[ht]))),gt++;if(pt>gt*.4&&mt/Math.abs(rt)>.25||pt>.9*gt&&gt>o){x=Math.max(H,r+o),Wt=!0;break}}if(Wt&&_<2){g=y,b=x,_++;continue}break}const M=Math.max(m/p,d/f),w=Math.min(m/p,d/f);let P=y,C=y,v=1/0,k=1/0;for(let L=y;L<=x;L++){const G=(c[Math.max(0,L-2)]+c[Math.max(0,L-1)]+c[L]+c[Math.min(i-1,L+1)]+c[Math.min(i-1,L+2)])/5,X=Math.abs(G-w-.1*(M-w)),j=Math.abs(G-w-.9*(M-w));X<v&&(v=X,P=L),j<k&&(k=j,C=L)}if(P<C){const L=P;P=C,C=L}const A=Math.max(4,Math.abs(P-C)*s),S=Math.max(a,a+2*Math.trunc(A/Math.max(s,1e-6)));P+=S,C-=S;const F=Math.max(Math.abs(P-r),Math.abs(C-r),Math.max(o,Math.trunc(4/Math.max(s,1e-6)))),T=new Array(i).fill(0),R=new Array(i).fill(0),E=1.85,O=.5;for(let L=0;L<n.length;L++){const G=Math.trunc(n[L]*It+r);let X=5;Math.abs(G-r)>O*F&&(X=Math.abs(G-r)>2*O*F?12:7);const j=Math.max(y,G-X),$=Math.min(x-1,G+X);if($<r-E*F||j>r+E*F){for(let et=j;et<=$;et++)T[et]+=t[L],R[et]+=1;continue}for(let et=j;et<=$;et++){let nt=1;if(Math.abs(et-r)<E*F){const ct=(et-r)*s;if(Math.abs(et-r)<F*O)nt=ci(n[L]-ct,li,s,1);else{const at=(Math.abs(et-r)/Math.max(1e-6,F)-O)/Math.max(1e-6,E-O),rt=1*(1-at)+.01*at;nt=ci(n[L]-ct,li,s,rt)}}!(nt>0)||!Number.isFinite(nt)||(T[et]+=t[L]*nt,R[et]+=nt)}}const U=new Array(i).fill(0);let I=-1,D=-1;for(let L=Math.max(0,y-1);L<=Math.min(i-1,x+1);L++)R[L]>0?(U[L]=T[L]/R[L],I<0&&(I=L),D=L):U[L]=An;if(I<0||D<0)return null;const B=3*It;let V=U[I],Q=1;for(let L=I+1;L<r&&Q<B;L++)U[L]!==An&&(V+=U[L],Q++);V/=Math.max(1,Q);let z=U[D],Y=1;for(let L=D-1;L>r&&Y<B;L--)U[L]!==An&&(z+=U[L],Y++);z/=Math.max(1,Y);for(let L=I-1;L>=0;L--)U[L]=V;for(let L=D+1;L<i;L++)U[L]=z;const W=Math.max(Math.trunc(r-E*F),y+2),K=Math.min(Math.trunc(r+E*F),x-3),it=Math.max(1,Math.trunc(2/Math.max(s,1e-6))),J=Jo(U,y,x,W,K,it),ot=new Array(i).fill(0);let Z=J[Math.max(0,Math.min(i-1,y))]??J[0]??0;for(let L=y;L<=x;L++){const G=J[L]??Z;ot[L]=(J[Math.min(i-1,L+1)]??G)-Z,Z=G}return{esf:J,lsfFull:ot}}function il(n){const t=new Array(n.length).fill(0);if(n.length===0)return t;t[0]=n[0];for(let e=1;e<n.length;e++){const i=pr(n[e]-n[e-1]);t[e]=t[e-1]+i}return t}function pr(n){if(!Number.isFinite(n))return n;let t=(n+Math.PI)%(2*Math.PI);return t<0&&(t+=2*Math.PI),t-Math.PI}function rl(n,t,e=0){if(n.length===0)return[];const i=Number.isFinite(e)?e:0,r=n.map((o,a)=>{const l=t[a]??0,u=-2*Math.PI*i*l;return pr(o-u)});return il(r).map((o,a)=>{const l=t[a]??0,u=-2*Math.PI*i*l;return o+u})}function sl(n,t,e,i=.05,r=Number.POSITIVE_INFINITY){const s=Math.min(n.length,t.length);if(s<2)return null;let o=0;if(e)for(let f=1;f<s;f++){const y=e[f];Number.isFinite(y)&&y>o&&(o=y)}const a=e&&o>0?Math.max(1e-6,o*i):0;let l=0,u=0,h=0,c=0,m=0,d=0;for(let f=1;f<s;f++){const y=n[f],x=t[f];if(!Number.isFinite(y)||!Number.isFinite(x)||Math.abs(y)<=1e-12||y>r)continue;const g=e?e[f]:1;if(!Number.isFinite(g)||g<=a)continue;const b=e?g*g:1;l+=b,u+=b*y,h+=b*x,c+=b*y*y,m+=b*y*x,d++}if(d<4||l<=0)return null;const p=l*c-u*u;return Math.abs(p)<=1e-12?null:{slope:(l*m-u*h)/p,intercept:(h*c-u*m)/p,used:d,threshold:a}}function al(n,t,e=Number.POSITIVE_INFINITY){const i=[],r=[],s=Math.min(n.length,t.length);for(let o=1;o<s;o++)Number.isFinite(n[o])&&Number.isFinite(t[o])&&Math.abs(n[o])>1e-12&&n[o]<=e&&(i.push(n[o]),r.push(t[o]));return i.length<2?{slope:0,intercept:Number.isFinite(t[0])?t[0]:0,used:i.length}:{...dr(i,r),used:i.length}}function ol(n,t,e,i=.05,r=Number.POSITIVE_INFINITY,s=0){const o=Math.min(n.length,t.length);if(o<4)return null;let a=0;if(e)for(let b=1;b<o;b++){const _=e[b];Number.isFinite(_)&&_>a&&(a=_)}const l=e&&a>0?Math.max(1e-6,a*i):0,u=[];for(let b=1;b<o;b++){const _=n[b],M=t[b];if(!Number.isFinite(_)||!Number.isFinite(M)||Math.abs(_)<=1e-12||_>r)continue;const w=e?e[b]:1;!Number.isFinite(w)||w<=l||u.push({freq:_,phase:M,weight:e?w*w:1})}if(u.length<4)return null;const h=b=>{let _=0,M=0;for(const k of u){const A=k.phase+2*Math.PI*b*k.freq;_+=k.weight*Math.sin(A),M+=k.weight*Math.cos(A)}const w=Math.atan2(_,M);let P=0,C=0;const v=.65;for(const k of u){const A=k.phase+2*Math.PI*b*k.freq,S=Math.abs(pr(A-w)),F=S<=v?S*S:v*(2*S-v);P+=k.weight*F,C+=k.weight}return{score:C>0?P/C:Number.POSITIVE_INFINITY,intercept:w}},c=Number.isFinite(s)?s:0,m=Math.max(2,Math.min(8,Math.abs(c)>1e-6?4:2)),d=.02;let p=c,f=h(p);for(let b=c-m;b<=c+m+d*.5;b+=d){const _=h(b);_.score<f.score&&(f=_,p=b)}let y=p-d*2,x=p+d*2;for(let b=0;b<32;b++){const _=y+(x-y)/3,M=x-(x-y)/3,w=h(_).score,P=h(M).score;w<P?x=M:y=_}const g=(y+x)*.5;return f=h(g),{slope:-2*Math.PI*g,intercept:f.intercept,used:u.length,threshold:l}}function ll(n,t){if(Number.isFinite(t)&&t>0)return t;let e=0;for(const i of n)Number.isFinite(i)&&i>e&&(e=i);return e>0?e:Number.POSITIVE_INFINITY}function mr(n,t,e,i=Number.POSITIVE_INFINITY,r=0){const s=rl(n,t,r),o=ll(t,i),a=ol(t,n,e,.05,o,r),l=a?null:sl(t,s,e,.05,o),u=a||l?null:al(t,s,o),h=(a==null?void 0:a.slope)??(l==null?void 0:l.slope)??(u==null?void 0:u.slope)??0,c=(a==null?void 0:a.intercept)??(l==null?void 0:l.intercept)??(u==null?void 0:u.intercept)??0,m=s.map((g,b)=>g-(h*(t[b]??0)+c)),d=Number.isFinite(m[0])?m[0]:0,p=t.map(g=>h*g+c+d),f=m.map(g=>g-d),y=c+d,x=Number.isFinite(h)?-h/(2*Math.PI):null;return{raw:[...n],unwrapped:s,linear:p,residual:f,fit:{groupDelayPx:x===null?null:x-r,absoluteGroupDelayPx:x,referenceDelayPx:r,slopeRadPerCycle:Number.isFinite(h)?h:null,interceptRad:Number.isFinite(y)?y:null,fitPointCount:(a==null?void 0:a.used)??(l==null?void 0:l.used)??(u==null?void 0:u.used)??0,fitWeightThreshold:(a==null?void 0:a.threshold)??(l==null?void 0:l.threshold)??0,fitDomain:"cycles-per-pixel",fitMaxFreqCyclesPerPixel:o}}}function gr(n,t,e){const i=[],r=[],s=[],o=[];for(let a=0;a<e.length;a++){const l=e[a];i.push(Ie(l,t,n.raw)),r.push(Ie(l,t,n.unwrapped)),s.push(Ie(l,t,n.linear)),o.push(Ie(l,t,n.residual))}return{ptfRaw:i,ptfUnwrapped:r,ptfLinear:s,ptfResidual:o}}function Bs(n,t){const e=n.map((a,l)=>({dist:a,value:t[l]})).filter(a=>Number.isFinite(a.dist)&&Number.isFinite(a.value)).sort((a,l)=>a.dist-l.dist);if(e.length===0)return{dists:[],vals:[]};const i=Math.max(1,Math.min(16,Math.floor(e.length*.1)));let r=0,s=0;for(let a=0;a<i;a++)r+=e[a].value,s+=e[e.length-1-a].value;r/=i,s/=i;const o=r<=s?e:e.map(a=>({dist:-a.dist,value:a.value})).sort((a,l)=>a.dist-l.dist);return{dists:o.map(a=>a.dist),vals:o.map(a=>a.value)}}function cl(n,t,e){const i=un(e,"RAW green-site sampling"),r=t%2,s=n%2;return i==="RGGB"||i==="BGGR"?(r+s)%2!==0:i==="GBRG"||i==="GRBG"?(r+s)%2===0:!1}function ul(n,t){return n+t&1?2:1}function es(n,t){return(t&1)<<1|n&1}function bn(n,t,e){return n===void 0?0:typeof n=="number"?Number.isFinite(n)?n:0:Number.isFinite(n[es(t,e)])?n[es(t,e)]:0}function _t(n,t,e,i){return i!==void 0&&i!=="default"?ul(n,t)===i:cl(n,t,e)}function ui(n){return n.length===0?0:n.reduce((t,e)=>t+e,0)/n.length}function Os(n,t,e){const i=(e%t+t)%t,r=Math.floor(i),s=(r+1)%t,o=i-r,a=r<n.length?n[r]:0,l=s<n.length?n[s]:0;return a*(1-o)+l*o}function zs(n,t){const e=n.length;if(e===0)return 0;if(t<=0)return n[0];if(t>=e-1)return n[e-1];const i=Math.floor(t),r=Math.min(e-1,i+1),s=t-i;return n[i]*(1-s)+n[r]*s}function hl(n,t,e){const i=n.length,r=new Array(i).fill(0);for(let s=0;s<i;s++)r[s]=zs(n,s-e+t);return r}function hi(n,t,e){const i=n.length;if(i===0)return{peakPos:0,peakIdx:0,peakVal:0};const r=Math.max(0,Math.floor(t-e)),s=Math.min(i-1,Math.ceil(t+e));let o=Math.max(0,Math.min(i-1,Math.round(t))),a=-1/0;for(let u=r;u<=s;u++){const h=Math.abs(n[u]);h>a&&(a=h,o=u)}Number.isFinite(a)||(a=Math.abs(n[o]??0));let l=o;if(o>0&&o<i-1){const u=n[o]>=0?1:-1,h=u*n[o-1],c=u*n[o],m=u*n[o+1],d=h-2*c+m;if(Number.isFinite(d)&&Math.abs(d)>1e-9){const p=.5*(h-m)/d;Number.isFinite(p)&&Math.abs(p)<=1&&(l=o+p)}}return{peakPos:l,peakIdx:o,peakVal:Math.abs(zs(n,l))}}function dl(n,t,e,i,r,s,o,a,l){const u=Math.floor(i.x),h=Math.floor(i.y),c=Math.floor(i.w),m=Math.floor(i.h),d=[],p=(f,y)=>{if(f<0||y<0||f>=t||y>=e)return null;const x=r+f,g=s+y;return Math.max(0,n[y*t+f]-bn(l,x,g))};for(let f=0;f<m;f++){const y=[],x=h+f;for(let g=0;g<c;g++){const b=u+g,_=r+b,M=s+x;if(_t(_,M,o,a)){y.push(p(b,x)??0);continue}const w=[],P=p(b-1,x),C=p(b+1,x),v=p(b,x-1),k=p(b,x+1);if(P!==null&&_t(_-1,M,o,a)&&w.push(P),C!==null&&_t(_+1,M,o,a)&&w.push(C),v!==null&&_t(_,M-1,o,a)&&w.push(v),k!==null&&_t(_,M+1,o,a)&&w.push(k),w.length===0){const A=[],S=p(b-1,x-1),F=p(b+1,x-1),T=p(b-1,x+1),R=p(b+1,x+1);S!==null&&_t(_-1,M-1,o,a)&&A.push(S),F!==null&&_t(_+1,M-1,o,a)&&A.push(F),T!==null&&_t(_-1,M+1,o,a)&&A.push(T),R!==null&&_t(_+1,M+1,o,a)&&A.push(R),y.push(ui(A));continue}y.push(ui(w))}d.push(y)}return d}function Vs(n,t,e,i,r,s,o,a,l){const u=Math.floor(i.x),h=Math.floor(i.y),c=Math.floor(i.w),m=Math.floor(i.h),d=[];for(let p=0;p<m;p++){const f=h+p;for(let y=0;y<c;y++){const x=u+y,g=r+x,b=s+f;_t(g,b,o,a)&&d.push({x:g,y:b,value:Math.max(0,n[f*t+x]-bn(l,g,b))})}}return d}function Gs(n,t,e,i,r){const s=Math.floor(i.x),o=Math.floor(i.y),a=Math.floor(i.w),l=Math.floor(i.h),u=(r==null?void 0:r.globalX)??0,h=(r==null?void 0:r.globalY)??0,c=!!(r!=null&&r.isThreePlane)&&n.length>=t*e*3,m=r==null?void 0:r.threePlaneChannel,d=[];for(let p=0;p<l;p++){const f=o+p,y=f*t;for(let x=0;x<a;x++){const g=s+x;let b=0;if(!c)b=Math.max(0,n[y+g]-bn(r==null?void 0:r.blackLevel,u+g,h+f));else{const _=(y+g)*3;if(m!==void 0)b=n[_+m];else{const M=n[_],w=n[_+1],P=n[_+2];b=.2126*M+.7152*w+.0722*P}}d.push({x:u+g,y:h+f,value:b})}}return d}function fl(n){var s;const t=n.length,e=((s=n[0])==null?void 0:s.length)??0;let i=0,r=0;for(let o=1;o<t-1;o++)for(let a=1;a<e-1;a++)i+=Math.abs(n[o][a+1]-n[o][a-1]),r+=Math.abs(n[o+1][a]-n[o-1][a]);return{gx:i,gy:r}}function Xs(n,t,e,i,r,s,o){var m;const a=n.length,l=((m=n[0])==null?void 0:m.length)??0,u=(d,p,f,y,x,g,b)=>{const _=b?a:l,M=Math.max(0,p-3),w=Math.min(_,p+4);let P=0,C=0;for(let k=M;k<w;k++)P+=d[k],C+=k*d[k];if(P<=0)return null;const v=C/P;return b?{x:t+y*x,y:f+v*g,weight:P}:{x:f+v*g,y:e+y*x,weight:P}},h=(d,p,f,y,x,g,b)=>{const _=Math.max(3,Math.min(Math.max(3,Math.floor(f/3)),Math.max(4,Math.round(f*.12)))),M=d.map(F=>{let T=-1/0,R=-1;for(let E=0;E<F.length;E++)F[E]>T&&(T=F[E],R=E);return{peakValue:T,peakIndex:R}}),w=(p-1)*.5,P=M.map((F,T)=>({...F,index:T})).filter(F=>F.peakValue>1&&F.peakIndex>=0).sort((F,T)=>{const R=T.peakValue-F.peakValue;return Math.abs(R)>1e-6?R:Math.abs(F.index-w)-Math.abs(T.index-w)});if(P.length===0)return[];const C=P[0],v=new Array(p).fill(null),k=u(d[C.index],C.peakIndex,y,C.index,x,g,b);if(!k)return[];v[C.index]=k;const A=(F,T)=>{const R=d[F],E=M[F];if(!(E.peakValue>1)||E.peakIndex<0)return null;const O=Math.max(0,Math.floor(T-_)),U=Math.min(R.length,Math.ceil(T+_+1));let I=-1/0,D=-1;for(let z=O;z<U;z++)R[z]>I&&(I=R[z],D=z);if(D<0||!(I>1))return null;const B=Math.max(1e-6,E.peakValue),V=Math.abs(D-T)<=_,Q=I>=B*.25;return!V||!Q?null:u(R,D,y,F,x,g,b)};let S=k?b?(k.y-y)/g:(k.x-y)/g:C.peakIndex;for(let F=C.index+1;F<p;F++){const T=A(F,S);T&&(v[F]=T,S=b?(T.y-y)/g:(T.x-y)/g)}S=k?b?(k.y-y)/g:(k.x-y)/g:C.peakIndex;for(let F=C.index-1;F>=0;F--){const T=A(F,S);T&&(v[F]=T,S=b?(T.y-y)/g:(T.x-y)/g)}return v.filter(F=>!!F)};if(s){const d=n.map(p=>p.map((f,y)=>y===0?0:Math.abs(f-p[y-1])));return h(d,a,l,t,r,i,!1)}const c=[];for(let d=0;d<l;d++){const p=new Array(a).fill(0);for(let f=1;f<a;f++)p[f]=Math.abs(n[f][d]-n[f-1][d]);c.push(p)}return h(c,l,a,e,i,r,!0)}function be(n){if(n.length<2)return null;let t=0,e=0,i=0;for(const c of n)t+=c.weight,e+=c.x*c.weight,i+=c.y*c.weight;if(t<=0)return null;e/=t,i/=t;let r=0,s=0,o=0;for(const c of n){const m=c.x-e,d=c.y-i;r+=c.weight*m*m,s+=c.weight*d*d,o+=c.weight*m*d}r/=t,s/=t,o/=t;const a=.5*Math.atan2(2*o,r-s);let l=Math.cos(a),u=Math.sin(a);const h=Math.hypot(l,u);return!Number.isFinite(h)||h<=1e-9?null:(l/=h,u/=h,(l<0||Math.abs(l)<=1e-9&&u<0)&&(l=-l,u=-u),{pointX:e,pointY:i,dirX:l,dirY:u,orientation:Math.abs(l)>=Math.abs(u)?1:2})}function pl(n,t){if(n.length!==4||t.length!==4||n.some(i=>i.length!==4))return null;const e=n.map((i,r)=>[...i,t[r]]);for(let i=0;i<4;i++){let r=i,s=Math.abs(e[i][i]);for(let a=i+1;a<4;a++){const l=Math.abs(e[a][i]);l>s&&(s=l,r=a)}if(!(s>1e-12))return null;if(r!==i){const a=e[i];e[i]=e[r],e[r]=a}const o=e[i][i];for(let a=i;a<=4;a++)e[i][a]/=o;for(let a=0;a<4;a++){if(a===i)continue;const l=e[a][i];if(!(Math.abs(l)<=1e-12))for(let u=i;u<=4;u++)e[a][u]-=l*e[i][u]}}return[e[0][4],e[1][4],e[2][4],e[3][4]]}function ml(n){if(n.length<4)return 0;const t=[...n].sort((m,d)=>m.x-d.x),e=t[0].x,r=t[t.length-1].x-e;if(!(r>1e-6))return 0;const s=16,o=[];for(let m=0;m<s;m++){const d=Math.max(0,Math.floor((m-1.5)*t.length/s)),p=Math.min(t.length-1,Math.floor((m+2.5)*t.length/s));if(p<d)continue;let f=0,y=0,x=0;for(let g=d;g<=p;g++)f+=t[g].x,y+=t[g].y,x++;x>0&&o.push({x:f/x,y:y/x})}if(o.length<4)return 0;const a=[.05952381,0,-.03571429,-.04761905,-.03571429,0,.05952381],l=new Array(o.length).fill(0),u=3;for(let m=0;m<o.length;m++){let d=0;for(let p=-u;p<=u;p++){const f=m+p;f<0||f>=o.length||(d+=a[p+u]*o[f].y)}l[m]=d}let h=0,c=1/0;for(let m=0;m<l.length-1;m++){const d=l[m],p=l[m+1];if(d===0){const b=Math.abs(o[m].x);b<c&&(c=b,h=o[m].x);continue}if(d*p>=0)continue;const f=p-d,y=Math.abs(f)>1e-12?-d/f:.5,x=o[m].x+(o[m+1].x-o[m].x)*y,g=Math.abs(x);g<c&&(c=g,h=x)}return!Number.isFinite(h)||h<e+.3*r||h>e+.7*r?0:h}function gl(n){if(n.length<8)return null;const t=[...n].filter(p=>Number.isFinite(p.x)&&Number.isFinite(p.y)&&Number.isFinite(p.weight)).sort((p,f)=>p.x-f.x);if(t.length<8)return null;const i=[ml(t),0,.5*(t[Math.floor((t.length-1)*.5)].x+t[Math.ceil((t.length-1)*.5)].x)];let r=null;for(const p of i){if(!Number.isFinite(p))continue;let f=0,y=0;for(const x of t)x.x<=p?f++:y++;if(f>=4&&y>=4){r=p;break}}if(r===null)return null;const s=Array.from({length:4},()=>new Array(4).fill(0)),o=new Array(4).fill(0);for(const p of t){const f=p.x,y=p.y,x=Math.max(1e-6,p.weight),g=f<=r?[f*f,f,1,0]:[2*r*f-r*r,f,1,(f-r)*(f-r)];for(let b=0;b<4;b++){o[b]+=x*g[b]*y;for(let _=0;_<4;_++)s[b][_]+=x*g[b]*g[_]}}const a=pl(s,o);if(!a)return null;const[l,u,h,c]=a,m=u+2*(l-c)*r,d=h+(c-l)*r*r;return[l,u,h,c,m,d].every(p=>Number.isFinite(p))?{splitX:r,left:[l,u,h],right:[c,m,d]}:null}function yl(n,t,e){const[i,r,s]=e;if(Math.abs(i)<=1e-12){const b=1+r*r;return b>1e-12?[(n-r*(s-t))/b]:[n]}const o=2*i*i,a=3*i*r,l=1+2*i*s-2*i*t+r*r,u=r*s-t*r-n;if(Math.abs(o)<=1e-12)return[n];const h=a/o,c=l/o,m=u/o,d=(h*h-3*c)/9,p=(2*h*h*h-9*h*c+27*m)/54,f=p*p-d*d*d;if(f<0&&d>0){const b=Math.acos(Math.max(-1,Math.min(1,p/Math.sqrt(d*d*d)))),_=-2*Math.sqrt(d);return[_*Math.cos(b/3)-h/3,_*Math.cos((b+2*Math.PI)/3)-h/3,_*Math.cos((b-2*Math.PI)/3)-h/3]}const y=Math.sqrt(Math.max(0,f)),x=-Math.sign(p||1)*Math.cbrt(Math.abs(p)+y),g=Math.abs(x)<=1e-12?0:d/x;return[x+g-h/3]}function Ys(n,t){if(n.length<8)return null;const e=-t.dirY,i=t.dirX,r=n.map(o=>({x:(o.x-t.pointX)*t.dirX+(o.y-t.pointY)*t.dirY,y:(o.x-t.pointX)*e+(o.y-t.pointY)*i,weight:o.weight})),s=gl(r);return s?{...t,normalX:e,normalY:i,splitX:s.splitX,left:s.left,right:s.right}:null}function xl(n,t){const e=n.x-t.pointX,i=n.y-t.pointY,r=e*t.dirX+i*t.dirY,s=e*t.normalX+i*t.normalY,o=r<t.splitX?t.left:t.right,a=yl(r,s,o);let l=s,u=Number.POSITIVE_INFINITY;for(const h of a){if(!Number.isFinite(h))continue;const c=o[0]*h*h+o[1]*h+o[2],m=r-h,d=s-c,p=Math.hypot(m,d);Number.isFinite(p)&&p<u&&(u=p,l=(d>=0?1:-1)*p)}return Number.isFinite(u)?l:s}function Ws(n,t,e,i,r,s,o,a=Tt){if(!t||t.length<8||n.length===0)return null;const l=t.filter(F=>Number.isFinite(F.x)&&Number.isFinite(F.y)).map(F=>({x:F.x,y:F.y,weight:1}));if(l.length<8)return null;const u=be(l);if(!u)return null;const h=e.p2.x-e.p1.x,c=e.p2.y-e.p1.y,m=Math.hypot(h,c);if(!Number.isFinite(m)||m<=1e-6)return null;let d=u.dirX,p=u.dirY;d*h+p*c<0&&(d=-d,p=-p);const f={...u,dirX:d,dirY:p},y=Ys(l,f),x=h/m,g=c/m,b=-g,_=x,M=(e.p1.x+e.p2.x)*.5,w=(e.p1.y+e.p2.y)*.5,P=-f.dirY,C=f.dirX,v=Math.abs(x)>=Math.abs(g)?1:2,k=[],A=[];for(const F of n){const T=F.x-M,R=F.y-w,E=T*x+R*g;if(Math.abs(E)>i)continue;const O=T*b+R*_;if(Math.abs(O)>r)continue;const U=y?xl(F,y):(F.x-f.pointX)*P+(F.y-f.pointY)*C;Number.isFinite(U)&&(k.push(U),A.push(F.value))}if(k.length<8)return null;const S=o?s!=null&&s.forceLegacyModel?wn(k,A,v,a,r*2):_n(k,A,v,a):bi(k,A,Math.max(2,r*2),s==null?void 0:s.manualBinSize,v,s==null?void 0:s.preferAutoPerEdgeBin);return S?(S.quadraticProjectionUsed=!!y,S):null}function ns(n){if(n.length<2)return null;const t=n.filter(e=>Number.isFinite(e.x)&&Number.isFinite(e.y)).map(e=>({x:e.x,y:e.y,weight:1}));return t.length<2?null:Tn(t,be(t))}function Tn(n,t){if(!t||n.length<2)return null;let e=1/0,i=-1/0;for(const s of n){const o=(s.x-t.pointX)*t.dirX+(s.y-t.pointY)*t.dirY;e=Math.min(e,o),i=Math.max(i,o)}if(!Number.isFinite(e)||!Number.isFinite(i))return null;const r=Math.max(.5,(i-e)*.03);return{p1:{x:t.pointX+t.dirX*(e-r),y:t.pointY+t.dirY*(e-r)},p2:{x:t.pointX+t.dirX*(i+r),y:t.pointY+t.dirY*(i+r)}}}function bl(n,t,e,i,r,s,o,a){return[Tn(n,t),Tn(e,i),Tn(r,s),Tn(o,a)]}function _l(n,t,e){if(!n||n.length<8)return;const i=n.filter(g=>Number.isFinite(g.x)&&Number.isFinite(g.y)).map(g=>({x:g.x,y:g.y,weight:1}));if(i.length<8)return;const r=be(i);if(!r)return;const s=t.p2.x-t.p1.x,o=t.p2.y-t.p1.y,a=Math.hypot(s,o);if(!Number.isFinite(a)||a<=1e-6)return;let l=r.dirX,u=r.dirY;l*s+u*o<0&&(l=-l,u=-u);const h=Ys(i,{...r,dirX:l,dirY:u});if(!h)return;const c=n.map(g=>(g.x-h.pointX)*h.dirX+(g.y-h.pointY)*h.dirY).filter(g=>Number.isFinite(g)),m=(t.p1.x-h.pointX)*h.dirX+(t.p1.y-h.pointY)*h.dirY,d=(t.p2.x-h.pointX)*h.dirX+(t.p2.y-h.pointY)*h.dirY;if(Number.isFinite(m)&&c.push(m),Number.isFinite(d)&&c.push(d),c.length<2)return;const p=Math.min(...c),f=Math.max(...c);if(!Number.isFinite(p)||!Number.isFinite(f)||f-p<=1e-6)return;const y=Math.max(21,e),x=[];for(let g=0;g<y;g++){const b=y===1?.5:g/(y-1),_=p+(f-p)*b,M=_<h.splitX?h.left:h.right,w=M[0]*_*_+M[1]*_+M[2];x.push({x:h.pointX+_*h.dirX+w*h.normalX,y:h.pointY+_*h.dirY+w*h.normalY})}return x}function wl(n,t,e,i){const r=Us(n,t);if(!r||r.profilePositionsPx.length<2||r.lsf.length<8)return null;let s=0,o=Number.POSITIVE_INFINITY;for(let a=0;a<r.profilePositionsPx.length;a++){const l=Math.abs(r.profilePositionsPx[a]);l<o&&(o=l,s=a)}return{esf:r.esf,lsfFull:r.lsf,binSize:r.diagnostics.profileStepPx,orientation:e,zeroIndex:s,shortSidePx:i,fallbackUsed:!1,continuousMtfV2:!0,continuousProjectedDistsPx:[...n],continuousProjectedValues:[...t],continuousMtfV2Result:r,mtfEngine:"continuous-v2"}}function Hs(n,t,e,i,r,s){if(n==="continuous-v2"){const l=wl(t,e,i,r);if(l)return l}const a=s();return a&&(a.mtfEngine="legacy"),a}function bi(n,t,e,i,r,s=!1,o=!1,a=!1){return Hs(a?"legacy":"continuous-v2",n,t,r,e,()=>Ml(n,t,e,i,r,s,o))}function Ml(n,t,e,i,r,s=!1,o=!1){if(n.length===0||t.length!==n.length)return null;const a=Bs(n,t),l=a.dists,u=a.vals;if(l.length===0)return null;const h=()=>{const M=e/2;let w=0;for(const C of n)Math.abs(C)<=M&&w++;if(w<=0)return .125;const P=40*M/w;return Math.max(.01,Math.min(.125,P))},c=(M,w,P,C,v)=>{if(!(C>0)||!(v>0)||!(P>w))return!1;const k=Math.floor((P-w)/C);if(k<2)return!1;const A=Math.max(w,-v),S=Math.min(P,v);if(!(S>A))return!1;const F=Math.max(0,Math.floor((A-w)/C)),T=Math.min(k,Math.ceil((S-w)/C));if(T<=F)return!1;const R=new Array(T-F).fill(0),E=w+F*C,O=w+T*C;for(let U=0;U<M.length;U++){const I=M[U];if(I<E)continue;if(I>=O)break;const D=Math.floor((I-w)/C);D>=F&&D<T&&R[D-F]++}return R.every(U=>U>0)},m=()=>{const M=l[0],w=l[l.length-1],P=Math.max(0,e*.25),C=.125,v=.5,k=.001,A=Math.round((v-C)/k);for(let S=0;S<=A;S++){const F=Number((C+S*k).toFixed(3));if(c(l,M,w,F,P))return F}return v};let d=.125;i&&i>0?d=Math.max(.01,Math.min(.5,i)):s?d=m():d=h();const p=l[0],f=l[l.length-1],y=Math.floor((f-p)/d);if(y<2)return null;const g=(()=>{const M=new Array(y).fill(0),w=new Array(y).fill(0);for(let C=0;C<l.length;C++){const v=(l[C]-p)/d;if(Number.isFinite(v))if(o){const k=Math.floor(v),A=v-k,S=1-A,F=A;k>=0&&k<y&&(M[k]+=u[C]*S,w[k]+=S);const T=k+1;T>=0&&T<y&&(M[T]+=u[C]*F,w[T]+=F)}else{const k=Math.floor(v);k>=0&&k<y&&(M[k]+=u[C],w[k]++)}}let P=u[0];for(let C=0;C<y;C++)w[C]>0?(M[C]/=w[C],P=M[C]):M[C]=P;return M})(),b=Ko(g,d),_=Math.max(0,Math.min(y-1,-p/d-.5));return{esf:g,lsfFull:b,binSize:d,orientation:r,zeroIndex:_,shortSidePx:e,fallbackUsed:!0}}function _n(n,t,e,i=Tt){if(n.length===0||t.length!==n.length)return null;const r=Bs(n,t),s=r.dists.map((l,u)=>({dist:l,value:r.vals[u]})).filter(l=>Math.abs(l.dist)<i);if(s.length<8)return null;const o=s.map(l=>l.dist),a=s.map(l=>l.value);return Hs("continuous-v2",o,a,e,i*2,()=>wn(o,a,e,i))}function wn(n,t,e,i=Tt,r=i*2){if(n.length===0||t.length!==n.length)return null;const s=[],o=[];for(let a=0;a<n.length;a++){const l=n[a],u=t[a];!Number.isFinite(l)||!Number.isFinite(u)||Math.abs(l)>=i||(s.push(l),o.push(u))}return s.length<8?null:bi(s,o,Math.max(2,r),void 0,e,!0,!0,!0)}function xe(n,t,e,i,r,s,o,a,l=Tt){if(s<=0||o<=0)return null;const h=!(!!(a!=null&&a.isThreePlane)&&n.length>=t*e*3)&&((a==null?void 0:a.greenOnly)??!1),c=h?un(a==null?void 0:a.bayerPattern,"constrained RAW SFR sampling"):null,m=r.p2.x-r.p1.x,d=r.p2.y-r.p1.y,p=Math.hypot(m,d);if(!Number.isFinite(p)||p<=1e-6)return null;const f=m/p,y=d/p,x=-y,g=f,b=(r.p1.x+r.p2.x)*.5,_=(r.p1.y+r.p2.y)*.5,M=Math.abs(f)>=Math.abs(y)?1:2,w=h?Vs(n,t,e,i,0,0,c,a==null?void 0:a.greenPhase,a==null?void 0:a.blackLevel):Gs(n,t,e,i,{...a,globalX:0,globalY:0});if(w.length===0)return null;if(!(a!=null&&a.disableQuadraticProjection)){const v=Ws(w,a==null?void 0:a.quadraticFitPoints,r,s,o,a,!0,l);if(v)return v}const P=[],C=[];for(const v of w){const k=v.x-b,A=v.y-_,S=k*f+A*y;if(Math.abs(S)>s)continue;const F=k*x+A*g;Math.abs(F)>o||(P.push(F),C.push(v.value))}return P.length<8?null:a!=null&&a.forceLegacyModel?wn(P,C,M,l):_n(P,C,M,l)}function Sl(n,t,e=0){const i=[...n.lsfFull];if(i.length<3)return!1;const r=Math.max(n.binSize,1e-6),s=Number.isFinite(n.zeroIndex)?n.zeroIndex:i.length/2,o=Math.max(1,Math.round((n.shortSidePx??0)*.5/r));let{peakPos:a,peakIdx:l,peakVal:u}=hi(i,s,o);const h=u*.2;let c=0,m=i.length-1;for(let f=l;f>=0;f--)if(i[f]<h){c=f;break}for(let f=l;f<i.length;f++)if(i[f]<h){m=f;break}const d=m-c;if(t&&d>0){const f=d*4,y=[],x=[];if(e>0){const g=Math.max(0,l-f-e),b=Math.max(0,l-f);for(let w=g;w<b;w++)y.push(w),x.push(i[w]);const _=Math.min(i.length,l+f),M=Math.min(i.length,l+f+e);for(let w=_;w<M;w++)y.push(w),x.push(i[w])}else{for(let g=0;g<Math.max(0,l-f);g++)y.push(g),x.push(i[g]);for(let g=Math.min(i.length,l+f);g<i.length;g++)y.push(g),x.push(i[g])}if(y.length>2){const{slope:g,intercept:b}=dr(y,x);for(let _=0;_<i.length;_++)i[_]=i[_]-(g*_+b);({peakPos:a}=hi(i,s,o))}}return Math.abs(a-s)*r<=Math.max(1e-6,(n.shortSidePx??0)/6)}function Pl(n){const t=n.length;if(t<3)return!1;let e=0,i=-1/0;for(let o=0;o<t;o++){const a=Math.abs(n[o]);a>i&&(i=a,e=o)}const r=t/3,s=2*t/3;return e>=r&&e<=s}function js(n){return n.mtfEngine==="continuous-v2"||n.continuousMtfV2===!0}function vl(n){return js(n)}function xt(n,t,e){const i=Math.max(0,Math.floor(n.x)),r=Math.max(0,Math.floor(n.y)),s=Math.min(t,Math.ceil(n.x+n.w)),o=Math.min(e,Math.ceil(n.y+n.h)),a=s-i,l=o-r;return a<2||l<2?null:{x:i,y:r,w:a,h:l}}function yr(n,t,e,i){const r=[],s=n.x,o=n.y,a=n.x+n.w,l=n.y+n.h,u=n.x+n.w*.5,h=n.y+n.h*.5,c=[{x:s,y:o},{x:a,y:o},{x:a,y:l},{x:s,y:l},{x:u,y:o},{x:a,y:h},{x:u,y:l},{x:s,y:h},{x:u,y:h}];for(const m of c){const d=we(m,t);Number.isFinite(d.x)&&Number.isFinite(d.y)&&r.push(d)}return r.length===0?null:xt(Et(r,2),e,i)}function Et(n,t=0){let e=1/0,i=1/0,r=-1/0,s=-1/0;for(const o of n)e=Math.min(e,o.x),i=Math.min(i,o.y),r=Math.max(r,o.x),s=Math.max(s,o.y);return{x:e-t,y:i-t,w:r-e+t*2,h:s-i+t*2}}function is(n,t){let e=Math.atan2(t,n)*180/Math.PI;return e<0&&(e+=180),e}function de(n,t){const e=n.p2.x-n.p1.x,i=n.p2.y-n.p1.y,r=Math.hypot(e,i);if(!Number.isFinite(r)||r<=1e-6)return null;const s=-i/r,o=e/r;return[{x:n.p1.x+s*t,y:n.p1.y+o*t},{x:n.p2.x+s*t,y:n.p2.y+o*t},{x:n.p2.x-s*t,y:n.p2.y-o*t},{x:n.p1.x-s*t,y:n.p1.y-o*t}]}function Cl(n,t,e,i,r,s,o,a){if(s<=0||o<=0)return null;const u=!(!!(a!=null&&a.isThreePlane)&&n.length>=t*e*3)&&((a==null?void 0:a.greenOnly)??!1),h=u?un(a==null?void 0:a.bayerPattern,"constrained RAW edge sampling"):null,c=r.p2.x-r.p1.x,m=r.p2.y-r.p1.y,d=Math.hypot(c,m);if(!Number.isFinite(d)||d<=1e-6)return null;const p=c/d,f=m/d,y=-f,x=p,g=(r.p1.x+r.p2.x)*.5,b=(r.p1.y+r.p2.y)*.5,_=Math.abs(p)>=Math.abs(f)?1:2,M=u?Vs(n,t,e,i,0,0,h,a==null?void 0:a.greenPhase,a==null?void 0:a.blackLevel):Gs(n,t,e,i,{...a,globalX:0,globalY:0});if(M.length===0)return null;if(!(a!=null&&a.disableQuadraticProjection)){const C=Ws(M,a==null?void 0:a.quadraticFitPoints,r,s,o,a,!1);if(C)return C}const w=[],P=[];for(const C of M){const v=C.x-g,k=C.y-b,A=v*p+k*f;if(Math.abs(A)>s)continue;const S=v*y+k*x;Math.abs(S)>o||(w.push(S),P.push(C.value))}return w.length<8?null:bi(w,P,o*2,a==null?void 0:a.manualBinSize,_,a==null?void 0:a.preferAutoPerEdgeBin)}function Fl(n,t,e){const i=[...n],r=new Array(n.length).fill(0),s=[0,0,0];let o=-1,a=1,l=-1;for(let c=1;c<n.length-1;c++){let m=0;if(n[c]>1e-4){m=Math.atan2(e[c]*o,t[c]*o);let d=0;for(let p=-5;p<=5;p++)Math.abs(m+p*2*Math.PI-s[1])<Math.abs(m+d*2*Math.PI-s[1])&&(d=p);m+=d*2*Math.PI}c>3&&Math.abs(m-s[0])>Math.PI/2&&l<c-1&&n[c]<.5&&(a*=-1,l=c),i[c]*=a,o*=-1,s[0]=s[1],s[1]=m,s[2]=m}const u=[-.086,.343,.486,.343,-.086];for(let c=0;c<n.length-3;c++){let m=0;for(let d=-2;d<=2;d++)m+=i[Math.abs(c+d)]*u[d+2];r[c]=m}for(let c=0;c<n.length-3;c++)i[c]=r[c];const h=7;for(let c=0;c<3;c++){r.fill(0);for(let d=0;d<n.length-h;d++)if(d<h)r[d]=i[d];else{const p=Math.min(5,Math.floor((d-5)/3)),f=$o[p];let y=0;for(let x=-h;x<=h;x++)y+=i[d+x]*f[x+h];r[d]=y}for(let d=n.length-h-2;d<n.length;d++)r[d]=i[d];const m=Math.abs(r[0])>1e-9?r[0]:1;for(let d=0;d<n.length;d++)i[d]=r[d]/m}for(let c=0;c<n.length;c++)i[c]=Math.abs(i[c]);return i}function kl(n,t){const e=[[0,0,0],[0,0,0],[0,0,0]],i=[0,0,0];for(let o=0;o<n.length;o++){const a=n[o],l=-t+o,u=[1,a,a*a];for(let h=0;h<3;h++){i[h]+=u[h]*l;for(let c=0;c<3;c++)e[h][c]+=u[h]*u[c]}}const r=la(e);if(!r)return null;const s=ca(r,i);return[s[0],s[1],s[2]]}function Al(n,t){let e=0,i=1,r=0,s=!1,o=0;const a=Math.min(n.length,ue/16*2);for(let l=0;l<a&&!s;l++){const u=n[l];if(i>.5&&u<=.5){const h=-(u-i)*ue;Math.abs(h)>1e-9&&(r=-((.5-i-h*e)/h),o=l,s=!0)}i=u,e=l/ue}if(!s)return null;if(o>=5&&o<a-10){const l=Math.min(Math.max(2,o-9),9),u=kl(n.slice(o-l,o+l+1),l);if(u){const c=(u[0]+.5*u[1]+.25*u[2]+o)/ue;if(o>9)r=c;else{const d=(o-5)/ue/8;r=(1-d)*r+d*c}}}return r*It*t}function Tl(n,t,e){for(let i=1;i<t.length;i++){if(t[i-1]<e||t[i]>=e)continue;const r=t[i-1]-t[i];if(Math.abs(r)<=1e-12)return n[i];const s=(t[i-1]-e)/r;return n[i-1]+s*(n[i]-n[i-1])}return null}function Il(n,t){const i=new gn(4096),r=new Float32Array(4096);let s=0;for(let f=1;f<n.lsf.length;f++)Math.abs(n.lsf[f])>Math.abs(n.lsf[s])&&(s=f);for(let f=0;f<4096;f++)r[f]=Os(n.lsf,4096,f+s);i.transform(r);const o=Math.max(1e-12,Math.hypot(i._real[0],i._imag[0])),a=Math.max(1e-9,n.diagnostics.profileStepPx),l=Math.max(...t),u=Math.min(4096/2,Math.ceil(l*4096*a)+2),h=[],c=[],m=[];for(let f=0;f<=u;f++)h.push(f/(4096*a)),c.push(f===0?0:Math.atan2(i._imag[f],i._real[f])),m.push(Math.hypot(i._real[f],i._imag[f])/o);const d=mr(c,h,m,l,0);return{...gr(d,h,t),fit:d.fit}}function Nl(n,t){const e=[];for(const m of n){const d=m.continuousMtfV2Result;if(d){e.push(d);continue}if(!m.continuousProjectedDistsPx||!m.continuousProjectedValues)continue;const p=Us(m.continuousProjectedDistsPx,m.continuousProjectedValues);p&&e.push(p)}if(e.length===0)return null;const i=new Array(501),r=new Array(501).fill(0);for(let m=0;m<=500;m++){const d=m/250;i[m]=d;for(const p of e)r[m]+=Ie(d,p.frequencies,p.mtf);r[m]/=e.length}const s=1,o=e.length===1?e[0].mtf50:Tl(i,r,.5),a=Math.min(...e.map(m=>m.diagnostics.maxReliableFrequencyCyclesPerPixel)),l=Array.from(new Set(e.flatMap(m=>m.diagnostics.warnings))),u=o===null?null:o<=a+1e-9,h=e[0],c=Il(h,i);return{esf:h.esf,lsf:[],lsfCropped:h.lsf,mtf:r,ptf:c.ptfResidual,ptfRaw:c.ptfRaw,ptfUnwrapped:c.ptfUnwrapped,ptfLinear:c.ptfLinear,ptfResidual:c.ptfResidual,ptfPhaseFit:c.fit,freqs:i.map(m=>m*s),mtf50:o===null?null:o*s,calcRadius:e.reduce((m,d)=>m+d.diagnostics.halfWidthPx,0)/e.length,mtfEngine:"continuous-v2",maxReliableFrequency:a*s,mtf50Reliable:u,warnings:l}}function Rl(n,t){if(n.length===0)return null;const e=ue,i=ue/16*4,r=new gn(e),s=1,o=el(),a=new Float32Array(501);for(let S=0;S<=500;S++)a[S]=S/500*s*2;const l=new Array(i).fill(0).map((S,F)=>F/e*s*It),u=new Float32Array(i).fill(0),h=new Float32Array(i).fill(0);let c=0,m=[],d=[],p=[],f=[],y=[],x=[],g=[],b=null,_=0;for(const S of n){const F=S.mtfmapperOrderedDists&&S.mtfmapperOrderedVals&&S.mtfmapperOrderedDists.length===S.mtfmapperOrderedVals.length?nl(S.mtfmapperOrderedDists,S.mtfmapperOrderedVals,S.mtfmapperEffectiveMaxDot??Tt):null,T=(F==null?void 0:F.lsfFull)??S.lsfFull,R=(F==null?void 0:F.esf)??S.esf;if(T.length<e)continue;const E=new Float32Array(e);for(let I=0;I<e;I++)E[I]=T[I]??0;r.transform(E);const O=Math.max(1e-9,Math.abs(r._real[0])),U=new Array(i).fill(0);for(let I=1;I<i;I++)U[I]=Math.atan2(r._imag[I],r._real[I]);for(let I=0;I<i;I++)u[I]+=r._real[I]/O,h[I]+=r._imag[I]/O;if(c++,_+=S.shortSidePx*.5,m.length===0){m=[...T],d=[...R];const I=new Array(i).fill(0);I[0]=1;for(let z=1;z<i;z++)I[z]=Math.hypot(r._real[z]/O,r._imag[z]/O);const D=l.map(z=>z),B=(Number.isFinite(S.zeroIndex)?S.zeroIndex:0)*(S.binSize??fr),V=mr(U,D,I,Number.POSITIVE_INFINITY,B),Q=gr(V,l,a);f=Q.ptfRaw,y=Q.ptfUnwrapped,x=Q.ptfLinear,g=Q.ptfResidual,p=Q.ptfResidual,b=V.fit}}if(c===0)return null;const M=new Float32Array(i),w=new Float32Array(i),P=new Array(i).fill(0);P[0]=1;for(let S=0;S<i;S++)M[S]=u[S]/c,w[S]=h[S]/c,S>0&&(P[S]=Math.hypot(M[S],w[S]));const C=Fl(P,M,w),v=new Array(i).fill(0);for(let S=0;S<i;S++)v[S]=C[S]/o[S];const k=Array.from(a,S=>Ie(S,l,v)),A=Al(v,s);return{esf:d,lsf:[],lsfCropped:m,mtf:k,ptf:p,ptfRaw:f,ptfUnwrapped:y,ptfLinear:x,ptfResidual:g,ptfPhaseFit:b,freqs:Array.from(a),mtf50:A,calcRadius:_/c}}function Ll(n){return n?{...n,mtfEngine:"legacy",maxReliableFrequency:null,mtf50Reliable:null,warnings:["Legacy MTF is active; a validated reliable-frequency limit is not available."]}:null}function El(n,t,e,i=!1,r=0,s=!1){if(n.length===0)return null;const o=n;if(!o||o.length===0)return null;if(n=o,n.every(js))return Nl(n);if(n.every(C=>C.mtfmapperLike))return Ll(Rl(n));const a=4096,l=new gn(a),u=1,h=new Float32Array(501);for(let C=0;C<=500;C++)h[C]=C/500*u*2;const c=new Float32Array(501).fill(0);let m=0,d=[],p=[],f=0,y=[],x=[],g=[],b=[],_=[],M=null;for(const C of n){let v=[...C.lsfFull];const k=C.binSize,A=Number.isFinite(C.zeroIndex)?C.zeroIndex:v.length/2,S=Math.max(1,Math.round((C.shortSidePx??0)*.5/Math.max(k,1e-6)));let{peakPos:F,peakIdx:T,peakVal:R}=hi(v,A,S);const E=R*.2;let O=0,U=v.length-1;for(let G=T;G>=0;G--)if(v[G]<E){O=G;break}for(let G=T;G<v.length;G++)if(v[G]<E){U=G;break}const I=U-O;let D=!1;if(i&&I>0){const G=I*4,X=[],j=[];if(r>0){const $=Math.max(0,T-G-r),et=Math.max(0,T-G);for(let at=$;at<et;at++)X.push(at),j.push(v[at]);const nt=Math.min(v.length,T+G),ct=Math.min(v.length,T+G+r);for(let at=nt;at<ct;at++)X.push(at),j.push(v[at])}else{for(let $=0;$<Math.max(0,T-G);$++)X.push($),j.push(v[$]);for(let $=Math.min(v.length,T+G);$<v.length;$++)X.push($),j.push(v[$])}if(X.length>2){const{slope:$,intercept:et}=dr(X,j);for(let nt=0;nt<v.length;nt++)v[nt]=v[nt]-($*nt+et);({peakPos:F,peakIdx:T,peakVal:R}=hi(v,A,S)),D=!0}}let B=0,V=0;if(t>0)V=t,B=Math.round(t/k);else{const G=R*.2;let X=0,j=v.length-1;for(let ct=T;ct>=0;ct--)if(v[ct]<G){X=ct;break}for(let ct=T;ct<v.length;ct++)if(v[ct]<G){j=ct;break}const et=(j-X)*k;let nt=Math.max(2,et*8);V=nt,B=Math.round(nt/k)}f+=V;const Q=Math.max(0,Math.floor(A-B)),z=Math.min(v.length,Math.ceil(A+B)),Y=v.slice(Q,z);if(Y.length<8)continue;const W=new Float32Array(a).fill(0),K=new Array(Y.length).fill(0);for(let G=0;G<Y.length;G++){let X=1;s&&(X=.5*(1-Math.cos(2*Math.PI*G/(Y.length-1)))),K[G]=Y[G]*X}const it=Math.max(0,Math.min(Y.length-1,F-Q));for(let G=0;G<a;G++)W[G]=Os(K,a,G+it);l.transform(W);const J=[],ot=[],Z=[];for(let G=0;G<=a/2;G++){const X=l._real[G],j=l._imag[G],$=Math.sqrt(X*X+j*j);J.push($),ot.push(G/(a*k)*u),Z.push(Math.atan2(j,X))}const L=J[0];if(L>0){for(let G=0;G<=500;G++){const X=h[G],$=tl(X,k);c[G]+=Ie(X,ot,J)/L/$}if(m++,d.length===0){d=hl(Y,it,(Y.length-1)/2),p=D?Ul(v):C.esf;const G=J.map(et=>et/L),X=ot.map(et=>et),j=mr(Z,X,G,Number.POSITIVE_INFINITY,0),$=gr(j,ot,h);x=$.ptfRaw,g=$.ptfUnwrapped,b=$.ptfLinear,_=$.ptfResidual,y=$.ptfResidual,M=j.fit}}}if(m===0)return null;const w=Array.from(c).map(C=>C/m);let P=null;for(let C=0;C<w.length-1;C++)if(w[C]>=.5&&w[C+1]<.5){P=h[C]+(.5-w[C])*(h[C+1]-h[C])/(w[C+1]-w[C]);break}return{esf:p,lsf:[],lsfCropped:d,mtf:w,ptf:y,ptfRaw:x,ptfUnwrapped:g,ptfLinear:b,ptfResidual:_,ptfPhaseFit:M,freqs:Array.from(h),mtf50:P,calcRadius:f/m,mtfEngine:"legacy",maxReliableFrequency:null,mtf50Reliable:null,warnings:["Legacy MTF is active; a validated reliable-frequency limit is not available."]}}function Ul(n){const t=new Array(n.length).fill(0);let e=0;for(let i=0;i<n.length;i++)e+=n[i],t[i]=e;return t}function Ie(n,t,e){if(n<=t[0])return e[0];if(n>=t[t.length-1])return e[e.length-1];let i=0;for(;n>t[i+1];)i++;const r=(n-t[i])/(t[i+1]-t[i]);return e[i]+r*(e[i+1]-e[i])}function xr(n){return{...Ji,...n,gradientPercentiles:n!=null&&n.gradientPercentiles&&n.gradientPercentiles.length>0?n.gradientPercentiles:Ji.gradientPercentiles}}function Dl(n){return!n||n.length===0?void 0:[Number.isFinite(n[0])?n[0]:0,Number.isFinite(n[1])?n[1]:Number.isFinite(n[0])?n[0]:0,Number.isFinite(n[2])?n[2]:Number.isFinite(n[0])?n[0]:0,Number.isFinite(n[3])?n[3]:Number.isFinite(n[0])?n[0]:0]}function Bl(n,t){const e=n.width,i=n.height,r=n.data,s=un(n.bayerPattern,"RAW SFR detection"),o=Dl(n.blackLevels),a=new Float32Array(e*i),l=(_,M)=>_<0||M<0||_>=e||M>=i?null:Math.max(0,r[M*e+_]-bn(o,_,M));let u=1/0,h=-1/0;for(let _=0;_<i;_++){const M=_*e;for(let w=0;w<e;w++){const P=M+w;let C=0;if(_t(w,_,s,t))C=l(w,_)??0;else{const v=[],k=l(w-1,_),A=l(w+1,_),S=l(w,_-1),F=l(w,_+1);if(k!==null&&_t(w-1,_,s,t)&&v.push(k),A!==null&&_t(w+1,_,s,t)&&v.push(A),S!==null&&_t(w,_-1,s,t)&&v.push(S),F!==null&&_t(w,_+1,s,t)&&v.push(F),v.length>0)C=ui(v);else{const T=[],R=l(w-1,_-1),E=l(w+1,_-1),O=l(w-1,_+1),U=l(w+1,_+1);R!==null&&_t(w-1,_-1,s,t)&&T.push(R),E!==null&&_t(w+1,_-1,s,t)&&T.push(E),O!==null&&_t(w-1,_+1,s,t)&&T.push(O),U!==null&&_t(w+1,_+1,s,t)&&T.push(U),C=ui(T)}}a[P]=C,C<u&&(u=C),C>h&&(h=C)}}if(!Number.isFinite(u)||!Number.isFinite(h)||h<=u+1e-9)return new Uint8Array(e*i);const c=1024,m=new Uint32Array(c),d=h-u;for(let _=0;_<a.length;_++){const M=Math.max(0,Math.min(1,(a[_]-u)/d)),w=Math.min(c-1,Math.max(0,Math.floor(M*(c-1))));m[w]++}const p=a.length,f=_=>{const M=p*_;let w=0;for(let P=0;P<c;P++)if(w+=m[P],w>=M)return u+P/Math.max(1,c-1)*d;return h},y=f(.01),x=f(.99),g=Math.max(1e-9,x-y),b=new Uint8Array(e*i);for(let _=0;_<a.length;_++){const M=Math.max(0,Math.min(1,(a[_]-y)/g));b[_]=Math.round(M*255)}return b}function Ol(n,t,e){const i=new Float32Array(n.length),r=new Float32Array(n.length),s=new Float32Array(n.length);for(let o=1;o<e-1;o++)for(let a=1;a<t-1;a++){const l=o*t+a,u=n[(o-1)*t+(a-1)],h=n[(o-1)*t+a],c=n[(o-1)*t+(a+1)],m=n[o*t+(a-1)],d=n[o*t+(a+1)],p=n[(o+1)*t+(a-1)],f=n[(o+1)*t+a],y=n[(o+1)*t+(a+1)],x=-u-2*m-p+(c+2*d+y),g=-u-2*h-c+(p+2*f+y);i[l]=x,r[l]=g,s[l]=Math.hypot(x,g)}return{gx:i,gy:r,magnitude:s}}function zl(n,t){let e=0,i=0;for(let l=0;l<n.length;l++){const u=n[l];!Number.isFinite(u)||u<=1e-6||(e=Math.max(e,u),i++)}if(i===0||e<=1e-6)return[];const r=1024,s=new Uint32Array(r);for(let l=0;l<n.length;l++){const u=n[l];if(!Number.isFinite(u)||u<=1e-6)continue;const h=Math.max(0,Math.min(1,u/e)),c=Math.min(r-1,Math.floor(h*(r-1)));s[c]++}const o=t&&t.length>0?t:Ji.gradientPercentiles,a=[];for(const l of o){const u=i*l;let h=0;for(let c=0;c<r;c++)if(h+=s[c],h>=u){a.push(c/Math.max(1,r-1)*e);break}}return Array.from(new Set(a.filter(l=>l>0))).sort((l,u)=>u-l)}function Vl(n,t){const e=new Uint8Array(n.length);for(let i=0;i<n.length;i++)e[i]=n[i]>=t?1:0;return e}const Gl=256*256;function Xl(n,t,e){if(n.length>=Gl){const s=Co.compute(n,t,e);if(s)return{gray:s.blurredGray,gradient:{gx:s.gx,gy:s.gy,magnitude:s.magnitude},backend:"webgl"}}const i=Zl(n,t,e),r=Ol(i,t,e);return{gray:i,gradient:r,backend:"cpu"}}function Yl(n,t,e,i){let r=n;for(let s=0;s<i;s++){const o=new Uint8Array(n.length);for(let a=0;a<e;a++)for(let l=0;l<t;l++){let u=0;for(let h=-1;h<=1&&!u;h++){const c=a+h;if(!(c<0||c>=e))for(let m=-1;m<=1;m++){const d=l+m;if(!(d<0||d>=t)&&r[c*t+d]){u=1;break}}}o[a*t+l]=u}r=o}return r}function Wl(n,t,e){const i=new Int32Array(n.length),r=[];let s=1;for(let o=0;o<n.length;o++){if(!n[o]||i[o]!==0)continue;const a=[o];i[o]=s;let l=0,u=t,h=e,c=0,m=0,d=0,p=!1;for(;l<a.length;){const f=a[l++],y=f%t,x=Math.floor(f/t);d++,u=Math.min(u,y),h=Math.min(h,x),c=Math.max(c,y),m=Math.max(m,x),(y===0||x===0||y===t-1||x===e-1)&&(p=!0);for(let g=-1;g<=1;g++)for(let b=-1;b<=1;b++){if(b===0&&g===0)continue;const _=y+b,M=x+g;if(_<0||M<0||_>=t||M>=e)continue;const w=M*t+_;!n[w]||i[w]!==0||(i[w]=s,a.push(w))}}r.push({label:s,x:u,y:h,w:c-u+1,h:m-h+1,area:d,touchesBorder:p}),s++}return{labels:i,components:r}}function qs(n,t){const e=Math.hypot(n,t);if(!Number.isFinite(e)||e<=1e-9)return null;let i=n/e,r=t/e;return(i<0||Math.abs(i)<=1e-9&&r<0)&&(i=-i,r=-r),{x:i,y:r}}function ve(n,t){if(n.length===0)return 0;const e=[...n].sort((o,a)=>o.value-a.value),i=e.reduce((o,a)=>o+Math.max(0,a.weight),0);if(i<=0)return e[Math.floor((e.length-1)*t)].value;const r=Math.max(0,Math.min(1,t))*i;let s=0;for(const o of e)if(s+=Math.max(0,o.weight),s>=r)return o.value;return e[e.length-1].value}function rs(n){const t=n.filter(i=>Number.isFinite(i)).sort((i,r)=>i-r);if(t.length===0)return 0;const e=i=>{if(i.length===1)return i[0];if(i.length===2)return(i[0]+i[1])*.5;const r=Math.ceil(i.length*.5);let s=0,o=1/0;for(let a=0;a+r-1<i.length;a++){const l=i[a+r-1]-i[a];l<o&&(o=l,s=a)}return e(i.slice(s,s+r))};return e(t)}function Hl(n,t,e,i,r,s,o){const a=[];for(let l=r.y;l<r.y+r.h;l++)for(let u=r.x;u<r.x+r.w;u++){const h=l*s+u;if(n[h]!==t||!e[h])continue;const c=i.magnitude[h];!Number.isFinite(c)||c<=1e-6||a.push({x:u,y:l,weight:c,gx:i.gx[h],gy:i.gy[h]})}return a}function jl(n){let t=0,e=0,i=0,r=0,s=0;for(const l of n){t+=l.weight,e+=l.x*l.weight,i+=l.y*l.weight;const u=Math.hypot(l.gx,l.gy);if(!Number.isFinite(u)||u<=1e-6)continue;const h=-l.gy/u,c=l.gx/u;r+=l.weight*(h*h-c*c),s+=l.weight*(2*h*c)}if(t<=0)return null;e/=t,i/=t;const o=.5*Math.atan2(s,r),a=qs(Math.cos(o),Math.sin(o));return a?{centerX:e,centerY:i,dirX:a.x,dirY:a.y,orthoX:-a.y,orthoY:a.x}:null}function Kn(n,t){let e=0,i=0;const r=-t.dirY,s=t.dirX;for(const o of n){const a=(o.x-t.pointX)*r+(o.y-t.pointY)*s;i+=o.weight*a*a,e+=o.weight}return e<=0?1/0:Math.sqrt(i/e)}function ss(n,t,e,i,r){const s=Math.max(0,Math.min(t-1,i)),o=Math.max(0,Math.min(e-1,r)),a=Math.floor(s),l=Math.floor(o),u=Math.min(t-1,a+1),h=Math.min(e-1,l+1),c=s-a,m=o-l,d=n[l*t+a],p=n[l*t+u],f=n[h*t+a],y=n[h*t+u],x=d+(p-d)*c,g=f+(y-f)*c;return x+(g-x)*m}function ql(n,t,e,i,r,s,o,a){const l=ss(n,t,e,i-s*a,r-o*a);return ss(n,t,e,i+s*a,r+o*a)-l}function $n(n,t,e,i,r){const s=Math.max(1e-6,e-t);if(n.length===0||!Number.isFinite(s))return{points:[],coverageRatio:0,centerCoverageRatio:0};const o=Math.max(1.5,Math.min(4,s/18)),a=Math.max(1,Math.ceil(s/o)),l=new Map;for(const p of n){const f=i(p);if(!Number.isFinite(f)||f<t||f>e)continue;const y=Math.max(0,Math.min(a-1,Math.floor((f-t)/o))),x=p.weight/(1+Math.abs(r(p))),g=l.get(y);(!g||x>g.score)&&l.set(y,{point:p,score:x})}const u=Array.from(l.values()).sort((p,f)=>i(p.point)-i(f.point)).map(p=>p.point),h=Math.max(0,Math.floor(a*.3)),c=Math.max(h+1,Math.ceil(a*.7));let m=0;for(let p=h;p<c;p++)l.has(p)&&m++;const d=Math.max(1,c-h);return{points:u,coverageRatio:u.length/a,centerCoverageRatio:m/d}}function Jn(n,t){const e=n.dirX*t.dirY-n.dirY*t.dirX;if(!Number.isFinite(e)||Math.abs(e)<=1e-6)return null;const i=t.pointX-n.pointX,r=t.pointY-n.pointY,s=(i*t.dirY-r*t.dirX)/e;return{x:n.pointX+n.dirX*s,y:n.pointY+n.dirY*s}}function Ql(n){if(n.length<3)return 0;let t=0;for(let e=0;e<n.length;e++){const i=n[e],r=n[(e+1)%n.length];t+=i.x*r.y-r.x*i.y}return t*.5}function Kl(n,t,e,i,r,s,o,a,l){const u=Hl(i,r.label,s,o,r,t),h=u.map(N=>({x:N.x,y:N.y}));if(u.length<l.minEdgePoints)return{candidate:null,failureStage:"min_edge_points",pointsCount:u.length,strongEdgePoints:h};const c=jl(u);if(!c)return{candidate:null,failureStage:"dominant_axes",pointsCount:u.length,strongEdgePoints:h};const m=u.map(N=>{const st=N.x-c.centerX,lt=N.y-c.centerY;return{...N,u:st*c.dirX+lt*c.dirY,v:st*c.orthoX+lt*c.orthoY}}),d={x:c.centerX,y:c.centerY},p=ve(m.map(N=>({value:N.u,weight:N.weight})),l.extentQuantileLow),f=ve(m.map(N=>({value:N.u,weight:N.weight})),l.extentQuantileHigh),y=ve(m.map(N=>({value:N.v,weight:N.weight})),l.extentQuantileLow),x=ve(m.map(N=>({value:N.v,weight:N.weight})),l.extentQuantileHigh),g=Math.max(1e-6,Math.max(Math.abs(p),Math.abs(f))),b=Math.max(1e-6,Math.max(Math.abs(y),Math.abs(x))),_=72,M=360/_,w=Array.from({length:_},()=>[]),P=N=>{let st=N%360;return st<0&&(st+=360),st},C=(N,st)=>{const lt=Math.abs(P(N)-P(st));return Math.min(lt,360-lt)};m.forEach(N=>{const st=N.u/g,lt=N.v/b,vt=P(Math.atan2(lt,st)*180/Math.PI),Bt=Math.hypot(st,lt),At=Math.max(0,Math.min(_-1,Math.floor(vt/M)));w[At].push({point:N,angleDeg:vt,normRadius:Bt})});const v=w.map(N=>N.length>0?rs(N.map(st=>st.normRadius)):-1/0),k=(N,st)=>{let lt=-1,vt=-1/0;for(let ft=0;ft<w.length;ft++){if(w[ft].length===0)continue;const ae=(ft+.5)*M;if(C(ae,N)>45||st.some(sa=>C(ae,sa)<45))continue;const Pe=v[ft];Pe>vt&&(vt=Pe,lt=ft)}let Bt=lt>=0?(lt+.5)*M:N,At=lt>=0?w[lt]:m.map(ft=>{const se=ft.u/g,ae=ft.v/b;return{point:ft,angleDeg:P(Math.atan2(ae,se)*180/Math.PI),normRadius:Math.hypot(se,ae)}}).filter(ft=>C(ft.angleDeg,N)<=45&&!st.some(se=>C(ft.angleDeg,se)<45));if(At.length===0&&(At=m.map(ft=>{const se=ft.u/g,ae=ft.v/b;return{point:ft,angleDeg:P(Math.atan2(ae,se)*180/Math.PI),normRadius:Math.hypot(se,ae)}}).filter(ft=>C(ft.angleDeg,N)<=45),Bt=N),At.length===0)return{x:m[0].x,y:m[0].y,u:m[0].u,v:m[0].v,angleDeg:N};const Wn=lt>=0?v[lt]:rs(At.map(ft=>ft.normRadius));let ge=0,Cn=0,Fr=0,kr=0,Ar=0;for(const ft of At){const se=C(ft.angleDeg,N)/45,ae=Math.abs(ft.normRadius-Wn),Pe=Math.max(1e-6,ft.point.weight)/(1+se*2+ae*6);ge+=Pe,Cn+=ft.point.x*Pe,Fr+=ft.point.y*Pe,kr+=ft.point.u*Pe,Ar+=ft.point.v*Pe}return ge>0?{x:Cn/ge,y:Fr/ge,u:kr/ge,v:Ar/ge,angleDeg:Bt}:{x:At[0].point.x,y:At[0].point.y,u:At[0].point.u,v:At[0].point.v,angleDeg:At[0].angleDeg}},A=k(225,[]),S=k(315,[A.angleDeg]),F=k(45,[A.angleDeg,S.angleDeg]),T=k(135,[A.angleDeg,S.angleDeg,F.angleDeg]),R=[{x:A.x,y:A.y},{x:S.x,y:S.y},{x:F.x,y:F.y},{x:T.x,y:T.y}],E=f-p,O=x-y,U=Math.min(E,O),I=Math.max(E,O);if(!Number.isFinite(U)||U<l.minSpanPx||I/Math.max(1,U)>l.maxAspectRatio)return{candidate:null,failureStage:"span_aspect",pointsCount:u.length,minSpan:U,maxSpan:I,axisCentroid:d,axisExtremePoints:R,strongEdgePoints:h};const D=Math.max(l.bandMinPx,Math.min(l.bandMaxPx,U*l.bandScale)),B=Math.max(1,Math.min(3,D*.55)),V=Math.max(a,0),Q=void 0,z=void 0,Y=N=>N.map(st=>({x:st.x,y:st.y,weight:st.weight})),W=N=>N.map(st=>({x:st.x,y:st.y})),K=(N,st,lt)=>N.filter(vt=>{if(!Number.isFinite(vt.weight)||vt.weight<V)return!1;const Bt=ql(n,t,e,vt.x,vt.y,st,lt,B);return Number.isFinite(Bt)&&Bt>=l.minPointContrast}),it=p,J=f,ot=y,Z=x,L=l.minCoverageRatio,G=l.minCenterCoverageRatio,X=[],j=[],$=[],et=[],nt=[],ct=(N,st,lt,vt,Bt,At)=>(lt-N)*(At-st)-(vt-st)*(Bt-N),at=N=>N>1e-6?1:N<-1e-6?-1:0,rt=[{u:(A.u+S.u)*.5,v:(A.v+S.v)*.5},{u:(S.u+F.u)*.5,v:(S.v+F.v)*.5},{u:(F.u+T.u)*.5,v:(F.v+T.v)*.5},{u:(T.u+A.u)*.5,v:(T.v+A.v)*.5}],Jt=(N,st)=>{const lt=at(ct(A.u,A.v,F.u,F.v,N,st)),vt=at(ct(S.u,S.v,T.u,T.v,N,st));return`${lt},${vt}`},Wt=new Map;rt.forEach((N,st)=>{Wt.set(Jt(N.u,N.v),st)});for(const N of m){if(!Number.isFinite(N.u)||!Number.isFinite(N.v)){nt.push(N);continue}let lt=Wt.get(Jt(N.u,N.v))??-1;if(lt<0){let vt=1/0;for(let Bt=0;Bt<rt.length;Bt++){const At=rt[Bt],Wn=(N.u-At.u)/g,ge=(N.v-At.v)/b,Cn=Wn*Wn+ge*ge;Cn<vt&&(vt=Cn,lt=Bt)}}lt===0?X.push(N):lt===1?j.push(N):lt===2?$.push(N):lt===3?et.push(N):nt.push(N)}const H=[...X,...$],pt=[...j,...et],mt={dir:H.length,ortho:pt.length,unassigned:m.length-H.length-pt.length},gt=X.length>=l.minSidePoints?ve(X.map(N=>({value:N.v,weight:N.weight})),.5):y,ht=$.length>=l.minSidePoints?ve($.map(N=>({value:N.v,weight:N.weight})),.5):x,Zt=et.length>=l.minSidePoints?ve(et.map(N=>({value:N.u,weight:N.weight})),.5):p,Me=j.length>=l.minSidePoints?ve(j.map(N=>({value:N.u,weight:N.weight})),.5):f,Le=[{x:(A.x+S.x)*.5,y:(A.y+S.y)*.5},{x:(S.x+F.x)*.5,y:(S.y+F.y)*.5},{x:(F.x+T.x)*.5,y:(F.y+T.y)*.5},{x:(T.x+A.x)*.5,y:(T.y+A.y)*.5}],Rt=X.filter(N=>Math.abs(N.v-gt)<=D&&N.u>=it&&N.u<=J),Ee=$.filter(N=>Math.abs(N.v-ht)<=D&&N.u>=it&&N.u<=J),Ue=et.filter(N=>Math.abs(N.u-Zt)<=D&&N.v>=ot&&N.v<=Z),Mn=j.filter(N=>Math.abs(N.u-Me)<=D&&N.v>=ot&&N.v<=Z),Ln=[Rt.length,Mn.length,Ee.length,Ue.length],En=[W(Rt),W(Mn),W(Ee),W(Ue)],Se=K(Rt,-c.orthoX,-c.orthoY),De=K(Ee,c.orthoX,c.orthoY),Be=K(Ue,-c.dirX,-c.dirY),pe=K(Mn,c.dirX,c.dirY),Sn=[Se.length,pe.length,De.length,Be.length],Pn=[W(Se),W(pe),W(De),W(Be)],te=$n(Se,it,J,N=>N.u,N=>N.v-gt),Oe=$n(pe,ot,Z,N=>N.v,N=>N.u-Me),ze=$n(De,it,J,N=>N.u,N=>N.v-ht),Ve=$n(Be,ot,Z,N=>N.v,N=>N.u-Zt),Ge=(N,st)=>N.slice().sort((lt,vt)=>st(lt)-st(vt)),Un=Ge(Se,N=>N.u),Dn=Ge(pe,N=>N.v),Bn=Ge(De,N=>N.u),On=Ge(Be,N=>N.v),_i=[Un.length,Dn.length,Bn.length,On.length],Vt=[te.coverageRatio,Oe.coverageRatio,ze.coverageRatio,Ve.coverageRatio];te.centerCoverageRatio,Oe.centerCoverageRatio,ze.centerCoverageRatio,Ve.centerCoverageRatio;const Ut=[W(Un),W(Dn),W(Bn),W(On)],Ht={axisPointCounts:mt,sideBandPointCounts:Ln,sideContrastPointCounts:Sn,gradientThreshold:a,pointAxisMinDot:Q,pointAxisMargin:z,bandWidth:D,minPointContrast:l.minPointContrast,minCoverageRatio:L,minCenterCoverageRatio:G,axisCentroid:d,axisExtremePoints:R,axisSideCenters:Le,strongEdgePoints:h,axisDirPoints:W(H),axisOrthoPoints:W(pt),axisUnassignedPoints:W(nt),sideBandPoints:En,sideContrastPoints:Pn};if(Un.length<l.minSidePoints||Bn.length<l.minSidePoints||On.length<l.minSidePoints||Dn.length<l.minSidePoints)return{candidate:null,failureStage:"min_side_points",pointsCount:u.length,minSpan:U,maxSpan:I,sidePointCounts:_i,sideCoverageRatios:Vt,...Ht,sideFitPoints:Ut};if(te.coverageRatio<L||Oe.coverageRatio<L||ze.coverageRatio<L||Ve.coverageRatio<L||te.centerCoverageRatio<G||Oe.centerCoverageRatio<G||ze.centerCoverageRatio<G||Ve.centerCoverageRatio<G)return{candidate:null,failureStage:"side_coverage",pointsCount:u.length,minSpan:U,maxSpan:I,sidePointCounts:_i,sideCoverageRatios:Vt,...Ht,sideFitPoints:Ut};const wi=Un,Mi=Bn,Si=On,Pi=Dn,jt=_i,Xe=be(Y(wi)),Ye=be(Y(Mi)),We=be(Y(Si)),He=be(Y(Pi));if(!Xe||!Ye||!We||!He)return{candidate:null,failureStage:"fit_lines",pointsCount:u.length,minSpan:U,maxSpan:I,sidePointCounts:jt,sideCoverageRatios:Vt,...Ht,sideFitPoints:Ut};const me=bl(Y(wi),Xe,Y(Pi),He,Y(Mi),Ye,Y(Si),We),zn=l.minAxisDot,Vn=(N,st,lt)=>Math.abs(N.dirX*st+N.dirY*lt),Dt=[Vn(Xe,c.dirX,c.dirY),Vn(He,c.orthoX,c.orthoY),Vn(Ye,c.dirX,c.dirY),Vn(We,c.orthoX,c.orthoY)];if(Dt[0]<zn||Dt[1]<zn||Dt[2]<zn||Dt[3]<zn)return{candidate:null,failureStage:"axis_alignment",pointsCount:u.length,minSpan:U,maxSpan:I,sidePointCounts:jt,sideCoverageRatios:Vt,axisDots:Dt,...Ht,sideFitPoints:Ut,sideFitLines:me};const qt=Math.max(l.residualLimitFloor,D*l.residualLimitScale),dt=[Kn(Y(wi),Xe),Kn(Y(Mi),Ye),Kn(Y(Si),We),Kn(Y(Pi),He)],vi=[dt[0],dt[3],dt[1],dt[2]],Ci=Math.max(...dt),ee=Jn(Xe,We),ne=Jn(Xe,He),ie=Jn(Ye,He),re=Jn(Ye,We);if(!ee||!ne||!ie||!re)return{candidate:null,failureStage:"corners",pointsCount:u.length,minSpan:U,maxSpan:I,sidePointCounts:jt,sideCoverageRatios:Vt,axisDots:Dt,sideResiduals:[dt[0],dt[3],dt[1],dt[2]],residualLimit:qt,...Ht,sideFitPoints:Ut,sideFitLines:me};const vn=[ee,ne,ie,re],Qt=Math.abs(Ql(vn));if(!Number.isFinite(Qt)||Qt<l.minQuadArea)return{candidate:null,failureStage:"quad_area",pointsCount:u.length,minSpan:U,maxSpan:I,sidePointCounts:jt,sideCoverageRatios:Vt,axisDots:Dt,sideResiduals:[dt[0],dt[3],dt[1],dt[2]],residualLimit:qt,...Ht,sideFitPoints:Ut,sideFitLines:me,quadArea:Qt};const Fi=Math.hypot(ne.x-ee.x,ne.y-ee.y),ki=Math.hypot(ie.x-ne.x,ie.y-ne.y),Ai=Math.hypot(ie.x-re.x,ie.y-re.y),Ti=Math.hypot(re.x-ee.x,re.y-ee.y),je=[Fi,ki,Ai,Ti],Ii=Math.min(Fi,ki,Ai,Ti),ea=Math.max(Fi,ki,Ai,Ti);if(!Number.isFinite(Ii)||Ii<l.minSideLength||ea/Math.max(1,Ii)>l.maxAspectRatio)return{candidate:null,failureStage:"side_length",pointsCount:u.length,minSpan:U,maxSpan:I,sidePointCounts:jt,sideCoverageRatios:Vt,axisDots:Dt,sideResiduals:[dt[0],dt[3],dt[1],dt[2]],residualLimit:qt,...Ht,sideFitPoints:Ut,sideFitLines:me,quadArea:Qt,sideLengths:je};const qe=qs(ne.x-ee.x+(ie.x-re.x),ne.y-ee.y+(ie.y-re.y));if(!qe)return{candidate:null,failureStage:"corners",pointsCount:u.length,minSpan:U,maxSpan:I,sidePointCounts:jt,sideCoverageRatios:Vt,axisDots:Dt,sideResiduals:[dt[0],dt[3],dt[1],dt[2]],residualLimit:qt,...Ht,sideFitPoints:Ut,sideFitLines:me,quadArea:Qt,sideLengths:je};const Pr={x:-qe.y,y:qe.x},Ni=(ee.x+ne.x+ie.x+re.x)*.25,Ri=(ee.y+ne.y+ie.y+re.y)*.25,Gn=vn.map(N=>{const st=N.x-Ni,lt=N.y-Ri;return{u:st*qe.x+lt*qe.y,v:st*Pr.x+lt*Pr.y}}),Xn=(Math.max(...Gn.map(N=>N.u))-Math.min(...Gn.map(N=>N.u)))*.5,Yn=(Math.max(...Gn.map(N=>N.v))-Math.min(...Gn.map(N=>N.v)))*.5;if(!Number.isFinite(Xn)||!Number.isFinite(Yn)||Math.min(Xn,Yn)<6)return{candidate:null,failureStage:"box_size",pointsCount:u.length,minSpan:U,maxSpan:I,sidePointCounts:jt,sideCoverageRatios:Vt,axisDots:Dt,sideResiduals:[dt[0],dt[3],dt[1],dt[2]],residualLimit:qt,sideFitPoints:Ut,quadArea:Qt,sideLengths:je};const q=ic(n,t,e,vn,Xe,He,Ye,We,Ni,Ri,Xn,Yn,l.innerPurityStdScale,l.outerMeanSpreadLimit);if(!Number.isFinite(Ci)||Ci>qt)return{candidate:null,failureStage:"residual",pointsCount:u.length,minSpan:U,maxSpan:I,sidePointCounts:jt,sideCoverageRatios:Vt,axisDots:Dt,sideResiduals:vi,residualLimit:qt,...Ht,sideFitPoints:Ut,sideFitLines:me,quadArea:Qt,sideLengths:je,outerContrast:q.contrast,outerUniformityOk:q.ok,outerMeanSpread:q.meanSpread,outerMeanSpreadLimit:q.meanSpreadLimit,outerAvgStd:q.avgStd,outerAvgStdLimit:q.avgStdLimit,outerSideMeans:q.outerSideMeans,outerSideStds:q.outerSideStds,outerSideStdLimit:q.outerSideStdLimit,outerSideQuads:q.outerSideQuads,innerSideUniformityOk:q.innerSideOk,innerSideStds:q.innerSideStds,innerSideStdLimit:q.innerSideStdLimit,innerSideQuads:q.innerSideQuads};const vr=l.filterBlockPurity&&(!q.ok||!q.innerSideOk);if(vr||q.contrast<l.minOuterContrast)return{candidate:null,failureStage:vr?q.ok?"inner_roi_uniformity":"outer_uniformity":"outer_contrast",pointsCount:u.length,minSpan:U,maxSpan:I,sidePointCounts:jt,sideCoverageRatios:Vt,axisDots:Dt,sideResiduals:vi,residualLimit:qt,...Ht,sideFitPoints:Ut,sideFitLines:me,quadArea:Qt,sideLengths:je,outerContrast:q.contrast,outerUniformityOk:q.ok,outerMeanSpread:q.meanSpread,outerMeanSpreadLimit:q.meanSpreadLimit,outerAvgStd:q.avgStd,outerAvgStdLimit:q.avgStdLimit,outerSideMeans:q.outerSideMeans,outerSideStds:q.outerSideStds,outerSideStdLimit:q.outerSideStdLimit,outerSideQuads:q.outerSideQuads,innerSideUniformityOk:q.innerSideOk,innerSideStds:q.innerSideStds,innerSideStdLimit:q.innerSideStdLimit,innerSideQuads:q.innerSideQuads};const Cr=xt(Et(vn,1),t,e);if(!Cr)return{candidate:null,failureStage:"bbox",pointsCount:u.length,minSpan:U,maxSpan:I,sidePointCounts:jt,sideCoverageRatios:Vt,axisDots:Dt,sideResiduals:[dt[0],dt[3],dt[1],dt[2]],residualLimit:qt,...Ht,sideFitPoints:Ut,sideFitLines:me,quadArea:Qt,sideLengths:je,outerContrast:q.contrast,outerUniformityOk:q.ok,outerMeanSpread:q.meanSpread,outerMeanSpreadLimit:q.meanSpreadLimit,outerAvgStd:q.avgStd,outerAvgStdLimit:q.avgStdLimit,outerSideMeans:q.outerSideMeans,outerSideStds:q.outerSideStds,outerSideStdLimit:q.outerSideStdLimit,outerSideQuads:q.outerSideQuads,innerSideUniformityOk:q.innerSideOk,innerSideStds:q.innerSideStds,innerSideStdLimit:q.innerSideStdLimit,innerSideQuads:q.innerSideQuads};const na=1/(1+Ci/Math.max(1,qt)),ia=l.filterBlockPurity?q.score:1,ra=q.contrast*ia*na*Math.sqrt(Qt);return{candidate:{centerX:Ni,centerY:Ri,dirX:qe.x,dirY:qe.y,halfWidth:Xn,halfHeight:Yn,score:ra,bbox:Cr,corners:vn,sideFitPoints:Ut,outerSideMeans:q.outerSideMeans,outerSideQuads:q.outerSideQuads},failureStage:null,pointsCount:u.length,minSpan:U,maxSpan:I,sidePointCounts:jt,sideCoverageRatios:Vt,axisDots:Dt,sideResiduals:vi,residualLimit:qt,...Ht,sideFitPoints:Ut,sideFitLines:me,quadArea:Qt,sideLengths:je,outerContrast:q.contrast,outerUniformityOk:q.ok,outerMeanSpread:q.meanSpread,outerMeanSpreadLimit:q.meanSpreadLimit,outerAvgStd:q.avgStd,outerAvgStdLimit:q.avgStdLimit,outerSideMeans:q.outerSideMeans,outerSideStds:q.outerSideStds,outerSideStdLimit:q.outerSideStdLimit,outerSideQuads:q.outerSideQuads,innerSideUniformityOk:q.innerSideOk,innerSideStds:q.innerSideStds,innerSideStdLimit:q.innerSideStdLimit,innerSideQuads:q.innerSideQuads}}function $l(n,t,e,i,r,s,o,a,l){return Kl(n,t,e,i,r,s,o,a,l).candidate}function Jl(n,t,e,i,r,s){const o=xr(r),a=Math.max(i*8,i+128),l=tc(n,t,e,o.downsampleMaxSide);s==null||s("Detecting candidates: downsampling...",.02),s==null||s("Detecting candidates: edge stage...",.06);const u=Xl(l.gray,l.width,l.height),h=u.gray,c=u.gradient;s==null||s(`Detecting candidates: gradient (${u.backend==="webgl"?"WebGL1":"CPU"})...`,.1);const m=zl(c.magnitude,o.gradientPercentiles),d=l.width*l.height,p=Math.max(o.minComponentAreaPx,Math.round(d*o.minComponentAreaRatio)),f=Math.max(p+1,Math.round(d*o.maxComponentAreaRatio)),y=[],x=Math.max(1,m.reduce((A,S,F)=>A+(F<=1,2),0));let g=0;for(let A=0;A<m.length;A++){const S=m[A],F=Vl(c.magnitude,S),T=A<=1?[3,2]:[2,1];for(const R of T){const E=g/x;s==null||s(`Detecting candidates: threshold ${A+1}/${m.length}, dilate ${R}`,.12+.78*E);const O=Yl(F,l.width,l.height,R),{labels:U,components:I}=Wl(O,l.width,l.height);for(const D of I){if(D.touchesBorder||D.area<p||D.area>f)continue;const B=$l(h,l.width,l.height,U,D,F,c,S,o);if(!B)continue;const V=1/l.scale,Q=B.corners.map(z=>({x:z.x*V,y:z.y*V}));y.push({centerX:B.centerX*V,centerY:B.centerY*V,dirX:B.dirX,dirY:B.dirY,halfWidth:B.halfWidth*V,halfHeight:B.halfHeight*V,score:B.score,bbox:{x:B.bbox.x*V,y:B.bbox.y*V,w:B.bbox.w*V,h:B.bbox.h*V},corners:Q,sideFitPoints:B.sideFitPoints?[B.sideFitPoints[0].map(z=>({x:z.x*V,y:z.y*V})),B.sideFitPoints[1].map(z=>({x:z.x*V,y:z.y*V})),B.sideFitPoints[2].map(z=>({x:z.x*V,y:z.y*V})),B.sideFitPoints[3].map(z=>({x:z.x*V,y:z.y*V}))]:void 0,outerSideMeans:B.outerSideMeans,outerSideQuads:B.outerSideQuads?[B.outerSideQuads[0].map(z=>({x:z.x*V,y:z.y*V})),B.outerSideQuads[1].map(z=>({x:z.x*V,y:z.y*V})),B.outerSideQuads[2].map(z=>({x:z.x*V,y:z.y*V})),B.outerSideQuads[3].map(z=>({x:z.x*V,y:z.y*V}))]:void 0})}y.length>a&&(y.sort((D,B)=>B.score-D.score),y.length=a),g++}}console.log(`[SFR Auto Detect] Candidate pool before dedupe: ${y.length}`),s==null||s(`Detecting candidates: deduplicating (0/${Math.max(1,Math.min(y.length,Math.max(i*4,i+32)))})...`,.94),y.sort((A,S)=>S.score-A.score);const b=Math.max(i*4,i+32),_=y.length>b?y.slice(0,b):y,M=[];if(_.length<=256){console.log(`[SFR Auto Detect] Using simple dedupe for ${_.length} candidates`);for(let A=0;A<_.length;A++){const S=_[A];console.log(`[SFR Auto Detect] Simple dedupe candidate ${A+1}/${_.length}`,S.bbox);const F=_.length<=0?1:A/_.length;if(s==null||s(`Detecting candidates: deduplicating (${A}/${_.length})...`,.94+.05*Math.min(1,F)),!M.some(R=>{const E=Math.hypot(S.centerX-R.centerX,S.centerY-R.centerY),O=Math.max(Math.hypot(S.bbox.w,S.bbox.h),Math.hypot(R.bbox.w,R.bbox.h));return as(S.bbox,R.bbox)>.28||E<O*.18})&&(M.push(S),M.length>=i))break}return s==null||s("Detecting candidates: deduplicating...",1),M}const w=Math.max(32,Math.round(Math.sqrt(Math.max(1,t*e)/4096))),P=new Map,C=new Set,v=A=>Math.floor(A/w),k=(A,S)=>{if(!Number.isFinite(A.bbox.x)||!Number.isFinite(A.bbox.y)||!Number.isFinite(A.bbox.w)||!Number.isFinite(A.bbox.h)||A.bbox.w<=0||A.bbox.h<=0||A.bbox.w>t*4||A.bbox.h>e*4)return;const F=v(A.bbox.x),T=v(A.bbox.x+A.bbox.w),R=v(A.bbox.y),E=v(A.bbox.y+A.bbox.h);for(let O=R;O<=E;O++)for(let U=F;U<=T;U++){const I=`${U},${O}`,D=P.get(I);D?D.push(S):P.set(I,[S])}};for(let A=0;A<_.length;A++){const S=_[A];if(A===0||A%200===0){const U=_.length<=0?1:A/_.length;s==null||s(`Detecting candidates: deduplicating (${A}/${_.length})...`,.94+.05*Math.min(1,U))}C.clear();const F=v(S.bbox.x),T=v(S.bbox.x+S.bbox.w),R=v(S.bbox.y),E=v(S.bbox.y+S.bbox.h);let O=!1;for(let U=R-1;U<=E+1&&!O;U++)for(let I=F-1;I<=T+1&&!O;I++){const D=P.get(`${I},${U}`);if(D)for(const B of D){if(C.has(B))continue;C.add(B);const V=M[B];if(!V)continue;const Q=Math.hypot(S.centerX-V.centerX,S.centerY-V.centerY),z=Math.max(Math.hypot(S.bbox.w,S.bbox.h),Math.hypot(V.bbox.w,V.bbox.h));if(as(S.bbox,V.bbox)>.28||Q<z*.18){O=!0;break}}}if(!O){const U=M.length;if(M.push(S),k(S,U),M.length>=i)break}}return s==null||s("Detecting candidates: deduplicating...",1),M}function Zl(n,t,e){const i=new Uint8Array(n.length);for(let r=0;r<e;r++)for(let s=0;s<t;s++){let o=0,a=0;for(let l=-1;l<=1;l++){const u=r+l;if(!(u<0||u>=e))for(let h=-1;h<=1;h++){const c=s+h;c<0||c>=t||(o+=n[u*t+c],a++)}}i[r*t+s]=Math.round(o/Math.max(1,a))}return i}function tc(n,t,e,i){const r=Math.max(t,e);if(r<=i)return{gray:n,width:t,height:e,scale:1};const s=i/r,o=Math.max(1,Math.round(t*s)),a=Math.max(1,Math.round(e*s)),l=new Uint8Array(o*a);for(let u=0;u<a;u++){const h=Math.min(e-1,Math.floor(u/s));for(let c=0;c<o;c++){const m=Math.min(t-1,Math.floor(c/s));l[u*o+c]=n[h*t+m]}}return{gray:l,width:o,height:a,scale:s}}function as(n,t){const e=Math.max(n.x,t.x),i=Math.max(n.y,t.y),r=Math.min(n.x+n.w,t.x+t.w),s=Math.min(n.y+n.h,t.y+t.h),o=Math.max(0,r-e),a=Math.max(0,s-i),l=o*a;if(l<=0)return 0;const u=n.w*n.h+t.w*t.h-l;return u>0?l/u:0}function os(n){const t=n.length;if(t===0)return{count:0,mean:0,std:1/0};let e=0;for(const s of n)e+=s;const i=e/t;let r=0;for(const s of n){const o=s-i;r+=o*o}return r/=t,{count:t,mean:i,std:Math.sqrt(Math.max(0,r))}}function ec(n,t,e,i){return{p1:{x:n.x-t*i,y:n.y-e*i},p2:{x:n.x+t*i,y:n.y+e*i}}}function ls(n,t,e,i,r){return[{x:n.p1.x+t*i,y:n.p1.y+e*i},{x:n.p2.x+t*i,y:n.p2.y+e*i},{x:n.p2.x+t*r,y:n.p2.y+e*r},{x:n.p1.x+t*r,y:n.p1.y+e*r}]}function nc(n,t,e){let i=0;for(let r=0;r<4;r++){const s=e[r],o=e[(r+1)%4],a=(o.x-s.x)*(t-s.y)-(o.y-s.y)*(n-s.x);if(Math.abs(a)<=1e-6)continue;const l=a>0?1:-1;if(i===0)i=l;else if(i!==l)return!1}return!0}function cs(n,t,e,i){const r=xt(Et(i,1),t,e);if(!r)return[];const s=[];for(let o=r.y;o<r.y+r.h;o++)for(let a=r.x;a<r.x+r.w;a++)nc(a,o,i)&&s.push(n[o*t+a]);return s}function ic(n,t,e,i,r,s,o,a,l,u,h,c,m,d){const p=h*2,f=c*2,y=Math.hypot(i[1].x-i[0].x,i[1].y-i[0].y),x=Math.hypot(i[2].x-i[1].x,i[2].y-i[1].y),g=Math.hypot(i[2].x-i[3].x,i[2].y-i[3].y),b=Math.hypot(i[3].x-i[0].x,i[3].y-i[0].y),M=Math.max(...[y,x,g,b]),w=Math.max(2,Math.min(p,f)),P=Math.max(4,M*.25),C=Math.max(2,Math.min(12,w*.22)),v=Math.max(1,Math.min(C,Math.max(1,w*.5-1))),k=1,A=Math.max(8,Math.round(Math.min(P,C*3))),S=[[i[0],i[1],i[1],i[0]],[i[1],i[2],i[2],i[1]],[i[2],i[3],i[3],i[2]],[i[3],i[0],i[0],i[3]]],F=[[i[0],i[1],i[1],i[0]],[i[1],i[2],i[2],i[1]],[i[2],i[3],i[3],i[2]],[i[3],i[0],i[0],i[3]]],T=[],R=[],E=[{corners:[i[0],i[1]],seedLine:r,sideLength:y},{corners:[i[1],i[2]],seedLine:s,sideLength:x},{corners:[i[2],i[3]],seedLine:o,sideLength:g},{corners:[i[3],i[0]],seedLine:a,sideLength:b}];for(let X=0;X<E.length;X++){const j=E[X],$=Math.max(1,j.sideLength*.5-1),et=Math.max(1,Math.min($,P*.5)),nt={x:(j.corners[0].x+j.corners[1].x)*.5,y:(j.corners[0].y+j.corners[1].y)*.5},ct=ec(nt,j.seedLine.dirX,j.seedLine.dirY,et),at=le(n,t,e,ct.p1,ct.p2,et,Math.max(4,k+Math.max(C,v)+2)),rt=(at==null?void 0:at.line)||ct,Jt=rt.p2.x-rt.p1.x,Wt=rt.p2.y-rt.p1.y,H=Math.hypot(Jt,Wt);if(!Number.isFinite(H)||H<=1e-6)return{ok:!1,score:0,meanSpread:1/0,meanSpreadLimit:1/0,avgStd:1/0,avgStdLimit:1/0,contrast:0,outerMean:0,outerSideMeans:[0,0,0,0],outerSideStds:[1/0,1/0,1/0,1/0],outerSideStdLimit:1/0,outerSideQuads:S,innerSideOk:!1,innerSideStds:[1/0,1/0,1/0,1/0],innerSideStdLimit:1/0,innerSideQuads:F};let pt=-Wt/H,mt=Jt/H;const gt={x:(rt.p1.x+rt.p2.x)*.5,y:(rt.p1.y+rt.p2.y)*.5};(gt.x-l)*pt+(gt.y-u)*mt<0&&(pt=-pt,mt=-mt);const ht=ls(rt,pt,mt,k,k+C),Zt=ls(rt,pt,mt,-k,-(k+v));S[X]=ht,F[X]=Zt,T.push(cs(n,t,e,ht)),R.push(cs(n,t,e,Zt))}const O=T.map(os);if(O.some(X=>X.count<A||!Number.isFinite(X.std)))return{ok:!1,score:0,meanSpread:1/0,meanSpreadLimit:1/0,avgStd:1/0,avgStdLimit:1/0,contrast:0,outerMean:0,outerSideMeans:[0,0,0,0],outerSideStds:[1/0,1/0,1/0,1/0],outerSideStdLimit:1/0,outerSideQuads:S,innerSideOk:!1,innerSideStds:[1/0,1/0,1/0,1/0],innerSideStdLimit:1/0,innerSideQuads:F};const U=R.map(os);if(U.some(X=>X.count<A||!Number.isFinite(X.std)||!Number.isFinite(X.mean)))return{ok:!1,score:0,meanSpread:1/0,meanSpreadLimit:1/0,avgStd:1/0,avgStdLimit:1/0,contrast:0,outerMean:0,outerSideMeans:[0,0,0,0],outerSideStds:[1/0,1/0,1/0,1/0],outerSideStdLimit:1/0,outerSideQuads:S,innerSideOk:!1,innerSideStds:[1/0,1/0,1/0,1/0],innerSideStdLimit:1/0,innerSideQuads:F};const I=O.map(X=>X.mean),D=I.reduce((X,j)=>X+j,0)/I.length,B=U.reduce((X,j)=>X+j.mean,0)/U.length,V=Math.abs(B-D),Q=Math.max(...I)-Math.min(...I),z=O.reduce((X,j)=>X+j.std,0)/O.length,Y=Math.max(0,d),W=Math.max(6,Math.min(20,V*.45)),K=I,it=O.map(X=>X.std),J=Math.max(W,Math.min(30,W*m)),ot=U.map(X=>X.std),Z=ot.every(X=>X<=J),L=Q<=Y&&z<=W,G=1/(1+Q/Math.max(1,Y)+z/Math.max(1,W));return{ok:L,score:G,meanSpread:Q,meanSpreadLimit:Y,avgStd:z,avgStdLimit:W,contrast:V,outerMean:D,outerSideMeans:K,outerSideStds:it,outerSideStdLimit:W,outerSideQuads:S,innerSideOk:Z,innerSideStds:ot,innerSideStdLimit:J,innerSideQuads:F}}function le(n,t,e,i,r,s,o){const a=r.x-i.x,l=r.y-i.y,u=Math.hypot(a,l);if(!Number.isFinite(u)||u<=1e-6)return null;const h=a/u,c=l/u,m=-c,d=h,p=(i.x+r.x)*.5,f=(i.y+r.y)*.5,y=de({p1:i,p2:r},o+2);if(!xt(Et(y||[i,r],2),t,e))return null;const g=Math.max(8,Math.round(s*2)+1),b=Math.max(8,Math.round(o*2)+1),_=g>1?s*2/(g-1):0,M=b>1?o*2/(b-1):0,w=Array.from({length:g},()=>new Array(b).fill(0));for(let S=0;S<g;S++){const F=-s+_*S;for(let T=0;T<b;T++){const R=-o+M*T,E=p+F*h+R*m,O=f+F*c+R*d;w[S][T]=rc(n,t,e,E,O)}}const P=Xs(w,-o,-s,M,_,!0);if(P.length<8)return null;const C=P.map(S=>{const F=S.x,T=S.y;return{x:p+T*h+F*m,y:f+T*c+F*d,weight:S.weight}}),v=be(C);if(!v)return null;let k=v.dirX,A=v.dirY;return k*h+A*c<0&&(k=-k,A=-A),{line:{p1:{x:v.pointX-k*s,y:v.pointY-A*s},p2:{x:v.pointX+k*s,y:v.pointY+A*s}},fitPoints:C.map(S=>({x:S.x,y:S.y}))}}function rc(n,t,e,i,r){if(t<=0||e<=0||n.length!==t*e)return 0;const s=Math.max(0,Math.min(t-1,i)),o=Math.max(0,Math.min(e-1,r)),a=Math.floor(s),l=Math.floor(o),u=Math.min(t-1,a+1),h=Math.min(e-1,l+1),c=s-a,m=o-l,d=n[l*t+a],p=n[l*t+u],f=n[h*t+a],y=n[h*t+u],x=d*(1-c)+p*c,g=f*(1-c)+y*c;return x*(1-m)+g*m}function sc(n,t,e,i,r){if(i<=0||r<=0||i>=t-1||r>=e-1)return{gx:0,gy:0};const s=r*t+i;return{gx:(n[s+1]-n[s-1])*.5,gy:(n[s+t]-n[s-t])*.5}}function ac(n){if(n.length<20)return null;const t=n.map(P=>Math.max(0,P.weight));let e=0;for(const P of t)e=Math.max(e,P);if(!(e>0))return null;for(let P=0;P<t.length;P++)t[P]/=e;const i=P=>{let C=0,v=0,k=0;for(let W=0;W<n.length;W++){const K=P[W];K>0&&(C+=K,v+=n[W].x*K,k+=n[W].y*K)}if(!(C>0))return null;v/=C,k/=C;let A=0,S=0,F=0;for(let W=0;W<n.length;W++){const K=P[W];if(!(K>0))continue;const it=n[W].x-v,J=n[W].y-k;A+=K*it*it,S+=K*it*J,F+=K*J*J}A/=C,S/=C,F/=C;const T=A+F,R=A*F-S*S,E=-T,O=R,U=Math.max(0,E*E-4*O),I=-.5*(E+(E>=0?1:-1)*Math.sqrt(U)),D=Math.abs(I)>1e-12?I:0,B=Math.abs(I)>1e-12?O/I:T,V=Math.max(D,B);let Q=0,z=1;Math.abs(S)>1e-10?(Q=V-F,z=S):A>F&&(Q=1,z=0);const Y=Math.atan2(-Q,z);return{centroid:{x:v,y:k},angle:Y,totalWeight:C}},r=i(t);if(!r)return null;const s=Math.cos(r.angle),o=Math.sin(r.angle),a=new Array(2*16*8).fill(0);for(let P=0;P<n.length;P++){const C=n[P].x-r.centroid.x,v=n[P].y-r.centroid.y,k=C*s+v*o,A=Math.round(k*8+16*8);if(A>=3&&A<a.length-3)for(let S=-3;S<=3;S++)a[A+S]+=t[P]}let l=16*8;for(let P=-5*8+16*8;P<=5*8+16*8;P++)a[P]>a[l]&&(l=P);let u=l-1;for(;u>1&&a[u]>.05*a[l];)u--;let h=l+1;for(;h<a.length-1&&a[h]>.05*a[l];)h++;let c=Math.max(1,u-8);for(;c>1&&a[c]<=a[u];)c--;let m=Math.min(a.length-1,h+8);for(;m<a.length-1&&a[m]<=a[h];)m++;const d=a.slice();for(let P=1;P<d.length;P++)d[P]+=d[P-1];const p=d[d.length-1];if(!(p>0))return null;let f=0;for(let P=1;P<d.length;P++)Math.abs(d[P]-.1*p)<Math.abs(d[f]-.1*p)&&(f=P);let y=d.length-1;for(let P=d.length-2;P>0;P--)Math.abs(d[P]-.9*p)<Math.abs(d[y]-.9*p)&&(y=P);let x=f/8-16,g=y/8-16;const b=g-x;x-=b*.7,g+=b*.7,x=Math.max((c+u)/16-16,x),g=Math.min((m+h)/16-16,g);const _=t.slice();for(let P=0;P<n.length;P++){const C=n[P].x-r.centroid.x,v=n[P].y-r.centroid.y,k=C*s+v*o;_[P]=k>=x&&k<=g?t[P]**4*(1/(10+Math.abs(k))):0}const M=i(_);if(!M)return null;const w=[];for(let P=0;P<n.length;P++)_[P]>0&&w.push({x:n[P].x,y:n[P].y,weight:_[P]});return w.length<8?null:{centroid:M.centroid,angle:M.angle,keptSamples:w}}function oc(n,t,e,i,r,s=Tt){var C;const o=r.x-i.x,a=r.y-i.y,l=Math.hypot(o,a);if(!Number.isFinite(l)||l<=12)return null;const u=o/l,h=a/l,c=5,m=4*s+.5,d=(v,k,A,S,F)=>{const T={x:v.x-k*l*.5,y:v.y-A*l*.5},R={p1:T,p2:{x:T.x+k*l,y:T.y+A*l}},E=de(R,m+2),O=xt(Et(E??[R.p1,R.p2],3),t,e),U=[],I=new Map;if(!O)return{reduced:null,scanlines:I};for(let D=O.y;D<O.y+O.h;D++)for(let B=O.x;B<O.x+O.w;B++){const V=B,Q=D,z=V-T.x,Y=Q-T.y,W=z*k+Y*A;if(!(W>c&&W<l-c))continue;const K=V-v.x,it=Q-v.y,J=K*S+it*F;if(Math.abs(J)<12){const{gx:ot,gy:Z}=sc(n,t,e,B,D),L=ot*ot+Z*Z;L>0&&U.push({x:V,y:Q,weight:L})}if(Math.abs(J)<m){const ot=I.get(D);ot?(B<ot.start&&(ot.start=B),B>ot.end&&(ot.end=B)):I.set(D,{start:B,end:B})}}return{reduced:ac(U),scanlines:I}};let p={x:(i.x+r.x)*.5,y:(i.y+r.y)*.5},f=u,y=h,x=-y,g=f,b=d(p,f,y,x,g);if(!b.reduced)return null;p=b.reduced.centroid,x=Math.cos(b.reduced.angle),g=Math.sin(b.reduced.angle),f=-g,y=x,f*u+y*h<0&&(f=-f,y=-y,x=-x,g=-g);let _=d(p,f,y,x,g);if(!_.reduced)return null;const M=Math.hypot(_.reduced.centroid.x-p.x,_.reduced.centroid.y-p.y);p=_.reduced.centroid,x=Math.cos(_.reduced.angle),g=Math.sin(_.reduced.angle),f=-g,y=x,f*u+y*h<0&&(f=-f,y=-y,x=-x,g=-g);let w=_;if(M>1){const v=d(p,f,y,x,g);v.reduced&&(w=v,p=v.reduced.centroid,x=Math.cos(v.reduced.angle),g=Math.sin(v.reduced.angle),f=-g,y=x,f*u+y*h<0&&(f=-f,y=-y))}const P=(((C=w.reduced)==null?void 0:C.keptSamples)??[]).map(v=>({x:v.x,y:v.y}));return P.length<8?null:{line:{p1:{x:p.x-f*l*.5,y:p.y-y*l*.5},p2:{x:p.x+f*l*.5,y:p.y+y*l*.5}},fitPoints:P,correctedScanlines:w.scanlines}}function lc(n,t,e){var b;const i=n.length,r=((b=n[0])==null?void 0:b.length)??0;if(r===0||i===0)return 0;const s=Math.max(0,Math.min(r-1,t)),o=Math.max(0,Math.min(i-1,e)),a=Math.floor(s),l=Math.floor(o),u=Math.min(r-1,a+1),h=Math.min(i-1,l+1),c=s-a,m=o-l,d=n[l][a],p=n[l][u],f=n[h][a],y=n[h][u],x=d*(1-c)+p*c,g=f*(1-c)+y*c;return x*(1-m)+g*m}function cc(n,t,e,i,r,s,o){var U;const a=i.p2.x-i.p1.x,l=i.p2.y-i.p1.y,u=Math.hypot(a,l);if(!Number.isFinite(u)||u<=1e-6)return null;const h=a/u,c=l/u,m=-c,d=h,p=(i.p1.x+i.p2.x)*.5,f=(i.p1.y+i.p2.y)*.5,y=de(i,s+2),x=xt(Et(y||[i.p1,i.p2],3),t,e);if(!x)return null;const g=dl(n,t,e,x,0,0,un(o==null?void 0:o.bayerPattern,"RAW edge refinement"),o==null?void 0:o.greenPhase,o==null?void 0:o.blackLevel),b=g.length;if((((U=g[0])==null?void 0:U.length)??0)<6||b<6)return null;const M=Math.max(8,Math.round(r*2)+1),w=Math.max(8,Math.round(s*2)+1),P=M>1?r*2/(M-1):0,C=w>1?s*2/(w-1):0,v=Array.from({length:M},()=>new Array(w).fill(0));for(let I=0;I<M;I++){const D=-r+P*I;for(let B=0;B<w;B++){const V=-s+C*B,Q=p+D*h+V*m,z=f+D*c+V*d;v[I][B]=lc(g,Q-x.x,z-x.y)}}const{gx:k,gy:A}=fl(v),S=k>=A,F=Xs(v,-s,-r,C,P,S);if(F.length<8)return null;const T=F.map(I=>{const D=I.x,B=I.y;return{x:p+B*h+D*m,y:f+B*c+D*d,weight:I.weight}}),R=be(T);if(!R)return null;let E=R.dirX,O=R.dirY;return E*h+O*c<0&&(E=-E,O=-O),{line:{p1:{x:R.pointX-E*r,y:R.pointY-O*r},p2:{x:R.pointX+E*r,y:R.pointY+O*r}},fitPoints:T.map(I=>({x:I.x,y:I.y}))}}function _e(n){const t=Math.max(0,Math.min(1,n));return t<=.04045?t/12.92:Math.pow((t+.055)/1.055,2.4)}function br(n,t,e){if(t<=0||e<=0||n.length!==t*e)return new Uint8Array(Math.max(0,t*e));let i=1/0,r=-1/0;for(let p=0;p<n.length;p++){const f=n[p];Number.isFinite(f)&&(f<i&&(i=f),f>r&&(r=f))}if(!Number.isFinite(i)||!Number.isFinite(r)||r<=i+1e-9)return new Uint8Array(t*e);const s=1024,o=new Uint32Array(s),a=r-i;for(let p=0;p<n.length;p++){const f=Math.max(0,Math.min(1,(n[p]-i)/a)),y=Math.min(s-1,Math.max(0,Math.floor(f*(s-1))));o[y]++}const l=n.length,u=p=>{const f=l*p;let y=0;for(let x=0;x<s;x++)if(y+=o[x],y>=f)return i+x/Math.max(1,s-1)*a;return r},h=u(.01),c=u(.99),m=Math.max(1e-9,c-h),d=new Uint8Array(t*e);for(let p=0;p<n.length;p++){const f=Math.max(0,Math.min(1,(n[p]-h)/m));d[p]=Math.round(f*255)}return d}function Zn(n,t,e=0){const i=new Float32Array(n.width*n.height),r=n.data;for(let s=0,o=0;s<r.length;s+=4,o++)i[o]=Qs(r,s,t,e);return br(i,n.width,n.height)}function uc(n){return Number.isFinite(n)?Math.max(0,Math.min(65535,Number(n))):0}function Qs(n,t,e,i=0){let r=n[t]/255,s=n[t+1]/255,o=n[t+2]/255;e&&(r=_e(r),s=_e(s),o=_e(o));const a=.2126*r+.7152*s+.0722*o;return Math.max(0,a-uc(i)/65535)}function hc(n,t){const e=n.width,i=n.height,r=n.data;if(r.length<e*i*3)return new Uint8Array(e*i);const s=new Float32Array(e*i);for(let o=0;o<e*i;o++){const a=o*3;t!==void 0?s[o]=r[a+t]:s[o]=.2126*r[a]+.7152*r[a+1]+.0722*r[a+2]}return br(s,e,i)}function dc(n){const t=new Float32Array(n.width*n.height);for(let e=0;e<n.data.length;e++)t[e]=n.data[e];return br(t,n.width,n.height)}function us(n,t,e){const i=xt(t,n.width,n.height);if(!i)return null;const r=new Uint16Array(i.w*i.h*3),s=n.data;let o=0;for(let a=i.y;a<i.y+i.h;a++)for(let l=i.x;l<i.x+i.w;l++){const u=(a*n.width+l)*4;let h=s[u]/255,c=s[u+1]/255,m=s[u+2]/255;e&&(h=_e(h),c=_e(c),m=_e(m)),r[o++]=Math.max(0,Math.min(65535,Math.round(h*65535))),r[o++]=Math.max(0,Math.min(65535,Math.round(c*65535))),r[o++]=Math.max(0,Math.min(65535,Math.round(m*65535)))}return{data:r,width:i.w,height:i.h}}function hs(n,t,e,i=0){const r=xt(t,n.width,n.height);if(!r)return null;const s=new Uint16Array(r.w*r.h*3),o=n.data;let a=0;for(let l=r.y;l<r.y+r.h;l++)for(let u=r.x;u<r.x+r.w;u++){const h=(l*n.width+u)*4,c=Math.max(0,Math.min(65535,Math.round(Qs(o,h,e,i)*65535)));s[a++]=c,s[a++]=c,s[a++]=c}return{data:s,width:r.w,height:r.h}}function fc(n,t){const e=xt(t,n.width,n.height);if(!e)return null;const i=new Uint16Array(e.w*e.h),r=n.data;let s=0;for(let o=e.y;o<e.y+e.h;o++)for(let a=e.x;a<e.x+e.w;a++){const l=(o*n.width+a)*4;i[s++]=Math.max(0,Math.min(65535,Math.round((.2126*r[l]+.7152*r[l+1]+.0722*r[l+2])*257)))}return{data:i,width:e.w,height:e.h}}function pc(n,t){const e=xt(t,n.width,n.height);if(!e)return null;const i=new Uint16Array(e.w*e.h);let r=0;for(let s=e.y;s<e.y+e.h;s++){const o=s*n.width;for(let a=e.x;a<e.x+e.w;a++)i[r++]=n.data[o+a]}return{data:i,width:e.w,height:e.h}}function Xt(n,t,e){return{x:n.x*t,y:n.y*e}}function mc(n,t,e){return{p1:Xt(n.p1,t,e),p2:Xt(n.p2,t,e)}}function Zi(n,t){const e=t(n);return{x:Number.isFinite(e.x)?e.x:n.x,y:Number.isFinite(e.y)?e.y:n.y}}function gc(n,t){return n.map(e=>Zi(e,t))}function oe(n,t,e,i=0,r=0){if(!n||n.length<8)return;const s=n.map(o=>({x:o.x*t-i,y:o.y*e-r})).filter(o=>Number.isFinite(o.x)&&Number.isFinite(o.y));return s.length>=8?s:void 0}function yc(n,t){return{p1:Zi(n.p1,t),p2:Zi(n.p2,t)}}function hn(n,t,e){return{p1:{x:n.p1.x-t,y:n.p1.y-e},p2:{x:n.p2.x-t,y:n.p2.y-e}}}function xc(n,t,e,i){const r=Math.max(0,Math.min(n.width-1,t)),o=(Math.max(0,Math.min(n.height-1,e))*n.width+r)*4;let a=n.data[o]/255,l=n.data[o+1]/255,u=n.data[o+2]/255;return i&&(a=_e(a),l=_e(l),u=_e(u)),(.2126*a+.7152*l+.0722*u)*65535}function Ks(n){return n.kind==="u16-mono"}function nn(n){return n.width}function rn(n){return n.height}function di(n,t,e,i){if(Ks(n)){const r=Math.max(0,Math.min(n.width-1,t)),s=Math.max(0,Math.min(n.height-1,e));return n.data[s*n.width+r]}return xc(n,t,e,i)}function bc(n,t,e,i){if(Ks(n)&&n.coordinateSpace==="distorted-padded"){const r=Math.round(n.paddingOffsetX??0),s=Math.round(n.paddingOffsetY??0);return di(n,t+r,e+s,i)}return di(n,t,e,i)}function _c(n,t,e,i,r=3){const o=[...n,{x:(n[0].x+n[1].x+n[2].x+n[3].x)*.25,y:(n[0].y+n[1].y+n[2].y+n[3].y)*.25},{x:(n[0].x+n[1].x)*.5,y:(n[0].y+n[1].y)*.5},{x:(n[1].x+n[2].x)*.5,y:(n[1].y+n[2].y)*.5},{x:(n[2].x+n[3].x)*.5,y:(n[2].y+n[3].y)*.5},{x:(n[3].x+n[0].x)*.5,y:(n[3].y+n[0].y)*.5}].map(a=>we(a,t)).filter(a=>Number.isFinite(a.x)&&Number.isFinite(a.y));return o.length===0?null:xt(Et(o,r),e,i)}function _r(n,t,e,i){const r=new Map;for(let s=n.y;s<n.y+n.h;s++)for(let o=n.x;o<n.x+n.w;o++){const a=Nt({x:o,y:s},t);if(!Number.isFinite(a.x)||!Number.isFinite(a.y))continue;const l=Math.round(a.x),u=Math.round(a.y);if(l<0||u<0||l>=e||u>=i)continue;const h=r.get(u);h?(l<h.start&&(h.start=l),l>h.end&&(h.end=l)):r.set(u,{start:l,end:l})}return r}function wc(n,t,e,i,r,s,o){const a=new Map,l=t.p2.x-t.p1.x,u=t.p2.y-t.p1.y,h=Math.hypot(l,u);if(!Number.isFinite(h)||h<=1e-6)return a;const c=l/h,m=u/h,d=-m,p=c,f={x:(t.p1.x+t.p2.x)*.5,y:(t.p1.y+t.p2.y)*.5},y=Math.max(1,e+1),x=Math.max(1,i+1.5);for(let g=n.y;g<n.y+n.h;g++)for(let b=n.x;b<n.x+n.w;b++){const _=b+.5,M=g+.5,w=_-f.x,P=M-f.y,C=w*c+P*m;if(!Number.isFinite(C)||Math.abs(C)>y)continue;const v=w*d+P*p;if(!Number.isFinite(v)||Math.abs(v)>x)continue;const k=Nt({x:_,y:M},r);if(!Number.isFinite(k.x)||!Number.isFinite(k.y))continue;const A=Math.round(k.x),S=Math.round(k.y);if(A<0||S<0||A>=s||S>=o)continue;const F=a.get(S);F?(A<F.start&&(F.start=A),A>F.end&&(F.end=A)):a.set(S,{start:A,end:A})}return a}function $s(n,t,e,i,r,s){const o=new Map,a=t.p2.x-t.p1.x,l=t.p2.y-t.p1.y,u=Math.hypot(a,l);if(!Number.isFinite(u)||u<=1e-6)return o;const h=a/u,c=l/u,m=-c,d=h,p={x:(t.p1.x+t.p2.x)*.5,y:(t.p1.y+t.p2.y)*.5},f=Math.max(1,e+1),y=Math.max(1,i+1.5),x=xt(n,r,s);if(!x)return o;for(let g=x.y;g<x.y+x.h;g++)for(let b=x.x;b<x.x+x.w;b++){const _=b-p.x,M=g-p.y,w=_*h+M*c;if(!Number.isFinite(w)||Math.abs(w)>f)continue;const P=_*m+M*d;if(!Number.isFinite(P)||Math.abs(P)>y)continue;const C=o.get(g);C?(b<C.start&&(C.start=b),b>C.end&&(C.end=b)):o.set(g,{start:b,end:b})}return o}function Js(n,t,e,i){const r=new Map;for(const[s,o]of n)for(let a=o.start;a<=o.end;a++){const l=we({x:a,y:s},t);if(!Number.isFinite(l.x)||!Number.isFinite(l.y))continue;const u=Math.round(l.x),h=Math.round(l.y);if(u<0||h<0||u>=e||h>=i)continue;const c=r.get(h);c?(u<c.start&&(c.start=u),u>c.end&&(c.end=u)):r.set(h,{start:u,end:u})}return r}function wr(n){return Math.abs(n.k1)<1e-4&&Math.abs(n.k2)<1e-4}function Mc(n){return[{x:n.x,y:n.y},{x:n.x+n.w,y:n.y},{x:n.x+n.w,y:n.y+n.h},{x:n.x,y:n.y+n.h}]}function Yt(n,t,e,i,r){return we({x:i.x+n*t,y:i.y+n*e},r)}function Mr(n,t,e,i,r){const o=Yt(n,t,e,i,r),a=Yt(n+1e-4,t,e,i,r);return{x:(a.x-o.x)/1e-4,y:(a.y-o.y)/1e-4}}function Zs(n,t,e,i,r,s){let o=.01;const a=h=>{const c=Yt(h,t,e,i,s);return Math.hypot(c.x-r.x,c.y-r.y)},l=a(n),u=a(n+o);if(!Number.isFinite(l)||!Number.isFinite(u))return null;if(l>u){let h=n,c=n+o;for(let m=0;m<24;m++){o*=2;const d=h+o,p=a(d),f=a(c);if(!Number.isFinite(p)||!Number.isFinite(f))break;if(p>=f)return{a:h,b:d};h=c,c=d}}else{let h=n,c=n+o;for(let m=0;m<24;m++){o*=2;const d=c-o,p=a(d),f=a(h);if(!Number.isFinite(p)||!Number.isFinite(f))break;if(p>=f)return{a:d,b:c};c=h,h=d}}return{a:n-Math.max(.5,o),b:n+Math.max(.5,o)}}function Sc(n,t,e=33){const i=n.p2.x-n.p1.x,r=n.p2.y-n.p1.y,s=Math.hypot(i,r);if(!Number.isFinite(s)||s<=1e-6)return[we(n.p1,t),we(n.p2,t)];const o=i/s,a=r/s,l={x:(n.p1.x+n.p2.x)*.5,y:(n.p1.y+n.p2.y)*.5},u=s*.5,h=Math.max(9,e),c=[];for(let m=0;m<h;m++){const d=h===1?.5:m/(h-1),p=-u+d*(u*2);c.push(Yt(p,o,a,l,t))}return c}function Pc(n,t,e,i,r,s,o=1){const a=n.p2.x-n.p1.x,l=n.p2.y-n.p1.y,u=Math.hypot(a,l);if(!Number.isFinite(u)||u<=1e-6)return null;const h=a/u,c=l/u,m={x:(n.p1.x+n.p2.x)*.5,y:(n.p1.y+n.p2.y)*.5},d=Math.max(24,Math.round(e*2)+1),p=[];for(let f=0;f<d;f++){const y=d===1?.5:f/(d-1),x=-e+y*(e*2),g=Yt(x,h,c,m,t),b=Mr(x,h,c,m,t),_=Math.hypot(b.x,b.y);if(!Number.isFinite(_)||_<=1e-9)continue;const M=-b.y/_,w=b.x/_;p.push({x:g.x+M*(i+o),y:g.y+w*(i+o)},{x:g.x-M*(i+o),y:g.y-w*(i+o)})}if(p.length<2){const f={p1:we(n.p1,t),p2:we(n.p2,t)},y=de(f,i+o);return y?xt(Et(y,o),r,s):null}return xt(Et(p,o),r,s)}function en(n,t,e,i,r){return Nt({x:i.x+n*t,y:i.y+n*e},r)}function vc(n,t,e,i,r,s){let o=.01;const a=h=>{const c=en(h,t,e,i,s);return Math.hypot(c.x-r.x,c.y-r.y)},l=a(n),u=a(n+o);if(!Number.isFinite(l)||!Number.isFinite(u))return null;if(l>u){let h=n,c=n+o;for(let m=0;m<24;m++){o*=2;const d=h+o,p=a(d),f=a(c);if(!Number.isFinite(p)||!Number.isFinite(f))break;if(p>=f)return{a:h,b:d};h=c,c=d}}else{let h=n,c=n+o;for(let m=0;m<24;m++){o*=2;const d=c-o,p=a(d),f=a(h);if(!Number.isFinite(p)||!Number.isFinite(f))break;if(p>=f)return{a:d,b:c};c=h,h=d}}return{a:n-Math.max(.5,o),b:n+Math.max(.5,o)}}function Sr(n,t,e){const i=(t.x-n.x)*(t.y-e.y)-(t.x-e.x)*(t.y-n.y);if(Math.abs(i)<=1e-12)return .5*(n.x+e.x);const r=(t.x-n.x)*(t.x-n.x)*(t.y-e.y)-(t.x-e.x)*(t.x-e.x)*(t.y-n.y),s=t.x-.5*r/i;return Number.isFinite(s)?s:.5*(n.x+e.x)}function Cc(n,t,e,i,r){const o=en(n,t,e,i,r),a=en(n+1e-4,t,e,i,r);return{x:(a.x-o.x)/1e-4,y:(a.y-o.y)/1e-4}}function Fc(n,t,e,i,r,s,o=!1,a){const l=[],u=[],h=a?Tt*2:Tt,c=Math.max(1,Math.min(s,h)),m=i.p2.x-i.p1.x,d=i.p2.y-i.p1.y,p=Math.hypot(m,d);if(!Number.isFinite(p)||p<=1e-6)return null;const f=m/p,y=d/p,x={x:(i.p1.x+i.p2.x)*.5,y:(i.p1.y+i.p2.y)*.5},g={p1:Nt(i.p1,e),p2:Nt(i.p2,e)},b=g.p2.x-g.p1.x,_=g.p2.y-g.p1.y,M=Math.hypot(b,_);if(!Number.isFinite(M)||M<=1e-6)return null;const w=b/M,P=_/M,C=-P,v=w,k={x:(g.p1.x+g.p2.x)*.5,y:(g.p1.y+g.p2.y)*.5},A=xt(a||Et(de(i,s+2)??[i.p1,i.p2],2),n.width,n.height);if(!A)return null;const S=wc(A,i,r,c,e,nn(t),rn(t));if(S.size===0)return null;const F=!wr(e);for(const[I,D]of S)if(!(I<0||I>=rn(t)))for(let B=D.start;B<=D.end;B++){if(B<0||B>=nn(t))continue;const V={x:B,y:I};let Q,z;if(F){const Y=we(V,e);if(!Number.isFinite(Y.x)||!Number.isFinite(Y.y)||Math.round(Y.x)<0||Math.round(Y.x)>=n.width||Math.round(Y.y)<0||Math.round(Y.y)>=n.height)continue;const W=Y.x-x.x,K=Y.y-x.y,it=W*f+K*y;if(!Number.isFinite(it))continue;Q=it,z=W*-y+K*f;const J=vc(it,f,y,x,V,e);if(!J)continue;const ot=.5*(J.a+J.b),Z=en(J.a,f,y,x,e),L=en(ot,f,y,x,e),G=en(J.b,f,y,x,e),X=Sr({x:J.a,y:Math.hypot(Z.x-V.x,Z.y-V.y)},{x:ot,y:Math.hypot(L.x-V.x,L.y-V.y)},{x:J.b,y:Math.hypot(G.x-V.x,G.y-V.y)});if(!Number.isFinite(X))continue;Q=X;const j=Cc(X,f,y,x,e),$=Math.hypot(j.x,j.y);if(!Number.isFinite($)||$<=1e-9)continue;const et=j.x/$,ct=-(j.y/$),at=et,rt=en(X,f,y,x,e);z=(V.x-rt.x)*ct+(V.y-rt.y)*at}else{const Y=V.x-k.x,W=V.y-k.y;Q=Y*w+W*P,z=Y*C+W*v}!Number.isFinite(Q)||Math.abs(Q)>r||!Number.isFinite(z)||Math.abs(z)>c||(l.push(z),u.push(di(t,B,I,o)))}if(l.length<8)return null;const T=Nt(i.p1,e),R=Nt(i.p2,e),E=R.x-T.x,O=R.y-T.y,U=Math.abs(E)>=Math.abs(O)?1:2;return _n(l,u,U,h)}function kc(n,t,e,i,r,s,o=!1,a,l,u,h=!1){const c=[],m=[],d=a?Tt*2:Tt,p=Math.max(1,Math.min(s,d)),f=i.p2.x-i.p1.x,y=i.p2.y-i.p1.y,x=Math.hypot(f,y);if(!Number.isFinite(x)||x<=1e-6)return null;const g=f/x,b=y/x,_=-b,M=g,w={x:(i.p1.x+i.p2.x)*.5,y:(i.p1.y+i.p2.y)*.5},P=xt(a||Et(de(i,d*4+2)??[i.p1,i.p2],2),n.width,n.height);if(!P)return null;const C=l??(u?_r(xt(u,nn(t),rn(t))??u,e,n.width,n.height):$s(P,i,Math.max(1,r),p*4+.5,n.width,n.height));if(C.size===0)return null;const v=Js(C,e,nn(t),rn(t));if(v.size===0)return null;const k=!wr(e);for(const[S,F]of v)for(let T=F.start;T<=F.end;T++){const R={x:T,y:S},E=Nt(R,e);if(!Number.isFinite(E.x)||!Number.isFinite(E.y)||Math.round(E.x)<0||Math.round(E.x)>=n.width||Math.round(E.y)<0||Math.round(E.y)>=n.height)continue;const O=E.x-w.x,U=E.y-w.y,I=O*g+U*b;let D=O*_+U*M;if(k){const B=Zs(I,g,b,w,R,e);if(!B)continue;const V=.5*(B.a+B.b),Q=Yt(B.a,g,b,w,e),z=Yt(V,g,b,w,e),Y=Yt(B.b,g,b,w,e),W=Sr({x:B.a,y:Math.hypot(Q.x-R.x,Q.y-R.y)},{x:V,y:Math.hypot(z.x-R.x,z.y-R.y)},{x:B.b,y:Math.hypot(Y.x-R.x,Y.y-R.y)});if(!Number.isFinite(W))continue;const K=Mr(W,g,b,w,e),it=Math.hypot(K.x,K.y);if(!Number.isFinite(it)||it<=1e-9)continue;const J=K.x/it,Z=-(K.y/it),L=J,G=Yt(W,g,b,w,e);D=(R.x-G.x)*Z+(R.y-G.y)*L}!Number.isFinite(I)||Math.abs(I)>Math.max(1,r)||!Number.isFinite(D)||Math.abs(D)>p||(c.push(D),m.push(bc(t,T,S,o)))}if(c.length<8)return null;const A=Math.abs(f)>=Math.abs(y)?1:2;return h?wn(c,m,A,d):_n(c,m,A,d)}function Ac(n,t,e,i,r,s){const o=n.width,a=n.height,l=un(n.bayerPattern,"corrected RAW edge measurement"),u=s!=null&&s.correctedRect?Tt*2:Tt,h=Math.max(1,Math.min(r,u)),c=(s==null?void 0:s.restrictToStrip)??!0,m=e.p2.x-e.p1.x,d=e.p2.y-e.p1.y,p=Math.hypot(m,d);if(!Number.isFinite(p)||p<=1e-6)return null;const f=m/p,y=d/p,x=-y,g=f,b={x:(e.p1.x+e.p2.x)*.5,y:(e.p1.y+e.p2.y)*.5},_={p1:{x:b.x-f*Math.max(1,i),y:b.y-y*Math.max(1,i)},p2:{x:b.x+f*Math.max(1,i),y:b.y+y*Math.max(1,i)}},M=de(_,h+2),w=(s!=null&&s.fixedRawRect?xt(s.fixedRawRect,o,a):null)??(s!=null&&s.correctedRect?yr(s.correctedRect,t,o,a):null)??(M?_c(M,t,o,a,2):null);if(!w)return null;const P=[],C=[];for(let A=w.y;A<w.y+w.h;A++){const S=A*o;for(let F=w.x;F<w.x+w.w;F++){if(!_t(F,A,l,s==null?void 0:s.greenPhase))continue;const T=Nt({x:F,y:A},t);if(!Number.isFinite(T.x)||!Number.isFinite(T.y))continue;const R=T.x-b.x,E=T.y-b.y,O=R*f+E*y;if(!Number.isFinite(O)||c&&Math.abs(O)>Math.max(1,i))continue;const U=R*x+E*g;if(!Number.isFinite(U)||c&&Math.abs(U)>h)continue;P.push(U);let I;I=Math.max(0,n.data[S+F]-bn(s==null?void 0:s.blackLevel,F,A)),C.push(I)}}if(P.length<8)return null;const v=Math.abs(m)>=Math.abs(d)?1:2,k=Math.max(2,(s==null?void 0:s.shortSidePxOverride)??(c?h*2:Math.min(w.w,w.h)));return bi(P,C,k,s==null?void 0:s.manualBinSize,v,s==null?void 0:s.preferAutoPerEdgeBin,!1,!!(s!=null&&s.forceLegacyModel))}function ta(n,t,e,i,r,s=!1,o,a,l,u=!1){if(a&&l){const C=l.p2.x-l.p1.x,v=l.p2.y-l.p1.y,k=Math.hypot(C,v);if(Number.isFinite(k)&&k>1e-6)return Fc(a,n,t,l,Math.max(1,k*.5),r,s,o)}const h=[],c=[],m=Tt,d=e.p2.x-e.p1.x,p=e.p2.y-e.p1.y,f=Math.hypot(d,p);if(!Number.isFinite(f)||f<=1e-6)return null;const y=d/f,x=p/f,g={x:(e.p1.x+e.p2.x)*.5,y:(e.p1.y+e.p2.y)*.5},b=(o?xt(o,nn(n),rn(n)):null)??Pc(e,t,i+1,r+1,nn(n),rn(n),1);if(!b)return null;const _=_r(b,t,nn(n),rn(n));if(_.size===0)return null;const M=-x,w=y;for(const[C,v]of _)for(let k=v.start;k<=v.end;k++){const A={x:k,y:C};let S=(A.x-g.x)*y+(A.y-g.y)*x,F=(A.x-g.x)*M+(A.y-g.y)*w;!Number.isFinite(S)||Math.abs(S)>i+1||!Number.isFinite(F)||Math.abs(F)>=m||(h.push(F),c.push(di(n,k,C,s)))}if(h.length<8)return null;const P=Math.abs(d)>=Math.abs(p)?1:2;return u?wn(h,c,P,m):_n(h,c,P,m)}function Tc(n,t,e){var m;const i=e.sourceMode??(t.isThreePlane?"three-plane":"rggb-raw"),r=e.useQuadraticProjection!==!1,s=!!e.forceRenderedMeasurement,o=n.width,a=n.height,l=e.threePlaneChannel,u=xr(e.detectionTuning),h=e.monochromeBlackLevel??0;if(i==="rggb-raw"&&!s){if(!t||t.isThreePlane)return null;const d=Bl(t,e.greenPhase),p=o/Math.max(1,t.width),f=a/Math.max(1,t.height);return{sourceMode:i,detectionGray:d,detectionWidth:t.width,detectionHeight:t.height,detectToMeasureX:1,detectToMeasureY:1,detectToDisplayX:p,detectToDisplayY:f,measureToDisplayX:p,measureToDisplayY:f,detectPointToDisplay:y=>Xt(y,p,f),measurePointToDisplay:y=>Xt(y,p,f),displayPointToDetect:y=>Xt(y,1/Math.max(1e-9,p),1/Math.max(1e-9,f)),measureUsesDisplayLine:!1,measureWidth:t.width,measureHeight:t.height,refineLine:(y,x,g)=>cc(t.data,t.width,t.height,{p1:y,p2:x},g*.5,Math.max(4,g*.2),{greenPhase:e.greenPhase,bayerPattern:t.bayerPattern})||le(d,t.width,t.height,y,x,g*.5,Math.max(4,g*.2)),measureEdge:(y,x,g,b,_)=>xe(t.data,t.width,t.height,y,x,g,b,{greenOnly:!0,greenPhase:e.greenPhase,bayerPattern:t.bayerPattern,blackLevel:e.blackLevel??void 0,preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?oe(_==null?void 0:_.fitPoints,1,1):void 0})}}if(s){const d=!!e.distortionCurveApplied&&!!e.distortionModel,p=i==="rggb-raw"&&!!e.distortionCorrected&&!!e.distortionModel&&!t.isThreePlane,f=!!e.distortionCorrected&&!!e.distortionModel&&!!e.distortionOriginalSamplingPlane,y=!!e.distortionCorrected&&!!e.distortionSamplingPlane,x=n,g=Zn(x,!!e.sfrHasGamma,i==="unmix-bw"?h:0);return{sourceMode:i,detectionGray:g,detectionWidth:x.width,detectionHeight:x.height,detectToMeasureX:1,detectToMeasureY:1,detectToDisplayX:1,detectToDisplayY:1,measureToDisplayX:1,measureToDisplayY:1,detectPointToDisplay:b=>b,measurePointToDisplay:b=>b,displayPointToDetect:b=>b,measureUsesDisplayLine:!1,measureWidth:x.width,measureHeight:x.height,refineLine:(b,_,M)=>(f?oc(g,n.width,n.height,b,_,Tt):null)||le(g,n.width,n.height,b,_,M*.5,Math.max(4,M*.2)),measureEdge:(b,_,M,w,P)=>{const C=e.distortionModel?yr(b,e.distortionModel,t.width,t.height):null;if(p){const S={p1:Nt(_.p1,e.distortionModel),p2:Nt(_.p2,e.distortionModel)},F=Math.hypot(S.p2.x-S.p1.x,S.p2.y-S.p1.y),T=Math.max(2,F*.5*u.sampleHalfWidthRatio);return Qo(t,e.distortionModel,S,Math.max(1,F*.5),T,{greenPhase:e.greenPhase,blackLevel:e.blackLevel??void 0,correctedRect:b})}if(p)return Ac(t,e.distortionModel,_,M,w,{greenPhase:e.greenPhase,blackLevel:e.blackLevel??void 0,correctedRect:b,fixedRawRect:C,preferAutoPerEdgeBin:!0});if(f)return kc(n,e.distortionOriginalSamplingPlane,e.distortionModel,_,M,w,!!e.sfrHasGamma,b,(P==null?void 0:P.correctedScanlines)??null,C);if(d){const S={p1:Nt(_.p1,e.distortionModel),p2:Nt(_.p2,e.distortionModel)},F=Math.hypot(S.p2.x-S.p1.x,S.p2.y-S.p1.y);return ta(e.distortionOriginalSamplingPlane??e.distortionSamplingPlane??e.distortionSamplingImage??n,e.distortionModel,S,Math.max(1,F*.5),w,!!e.sfrHasGamma,b,e.distortionBaseImage??n,_)}if(y){const S=pc(e.distortionSamplingPlane,b);if(!S)return null;const F=hn(_,b.x,b.y);return xe(S.data,S.width,S.height,{x:0,y:0,w:S.width,h:S.height},F,M,w,{preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?oe(P==null?void 0:P.fitPoints,1,1,b.x,b.y):void 0})}const k=i==="unmix-bw"?hs(n,b,!!e.sfrHasGamma,h):us(n,b,!!e.sfrHasGamma);if(!k)return null;const A=hn(_,b.x,b.y);return xe(k.data,k.width,k.height,{x:0,y:0,w:k.width,h:k.height},A,M,w,{isThreePlane:!0,threePlaneChannel:void 0,preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?oe(P==null?void 0:P.fitPoints,1,1,b.x,b.y):void 0})}}}if(i==="three-plane"){if(t.isThreePlane&&!e.sfrHasGamma){const p=hc(t,l),f=o/Math.max(1,t.width),y=a/Math.max(1,t.height);return{sourceMode:i,detectionGray:p,detectionWidth:t.width,detectionHeight:t.height,detectToMeasureX:1,detectToMeasureY:1,detectToDisplayX:f,detectToDisplayY:y,measureToDisplayX:f,measureToDisplayY:y,detectPointToDisplay:x=>Xt(x,f,y),measurePointToDisplay:x=>Xt(x,f,y),displayPointToDetect:x=>Xt(x,1/Math.max(1e-9,f),1/Math.max(1e-9,y)),measureUsesDisplayLine:!1,measureWidth:t.width,measureHeight:t.height,refineLine:(x,g,b)=>le(p,t.width,t.height,x,g,b*.5,Math.max(4,b*.2)),measureEdge:(x,g,b,_,M)=>xe(t.data,t.width,t.height,x,g,b,_,{isThreePlane:!0,threePlaneChannel:l,preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?oe(M==null?void 0:M.fitPoints,1,1):void 0})}}const d=Zn(n,!!e.sfrHasGamma);return{sourceMode:i,detectionGray:d,detectionWidth:n.width,detectionHeight:n.height,detectToMeasureX:1,detectToMeasureY:1,detectToDisplayX:1,detectToDisplayY:1,measureToDisplayX:1,measureToDisplayY:1,detectPointToDisplay:p=>p,measurePointToDisplay:p=>p,displayPointToDetect:p=>p,measureUsesDisplayLine:!1,measureWidth:n.width,measureHeight:n.height,refineLine:(p,f,y)=>le(d,n.width,n.height,p,f,y*.5,Math.max(4,y*.2)),measureEdge:(p,f,y,x,g)=>{const b=us(n,p,!!e.sfrHasGamma);if(!b)return null;const _=hn(f,p.x,p.y);return xe(b.data,b.width,b.height,{x:0,y:0,w:b.width,h:b.height},_,y,x,{isThreePlane:!0,threePlaneChannel:l,preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?oe(g==null?void 0:g.fitPoints,1,1,p.x,p.y):void 0})}}}if(i==="unmix-bw"){if(t&&!t.isThreePlane&&e.displaySettings){const p=oo(t,e.displaySettings,e.blackLevel??e.monochromeBlackLevel??void 0);if(p){const f=dc(p),y=o/Math.max(1,t.width),x=a/Math.max(1,t.height);return{sourceMode:i,detectionGray:f,detectionWidth:t.width,detectionHeight:t.height,detectToMeasureX:1,detectToMeasureY:1,detectToDisplayX:y,detectToDisplayY:x,measureToDisplayX:y,measureToDisplayY:x,detectPointToDisplay:g=>Xt(g,y,x),measurePointToDisplay:g=>Xt(g,y,x),displayPointToDetect:g=>Xt(g,1/Math.max(1e-9,y),1/Math.max(1e-9,x)),measureUsesDisplayLine:!1,measureWidth:t.width,measureHeight:t.height,refineLine:(g,b,_)=>le(f,t.width,t.height,g,b,_*.5,Math.max(4,_*.2)),measureEdge:(g,b,_,M,w)=>xe(p.data,p.width,p.height,g,b,_,M,{preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?oe(w==null?void 0:w.fitPoints,1,1):void 0})}}}const d=Zn(n,!!e.sfrHasGamma,h);return{sourceMode:i,detectionGray:d,detectionWidth:n.width,detectionHeight:n.height,detectToMeasureX:1,detectToMeasureY:1,detectToDisplayX:1,detectToDisplayY:1,measureToDisplayX:1,measureToDisplayY:1,detectPointToDisplay:p=>p,measurePointToDisplay:p=>p,displayPointToDetect:p=>p,measureUsesDisplayLine:!1,measureWidth:n.width,measureHeight:n.height,refineLine:(p,f,y)=>le(d,n.width,n.height,p,f,y*.5,Math.max(4,y*.2)),measureEdge:(p,f,y,x,g)=>{const b=hs(n,p,!!e.sfrHasGamma,h);if(!b)return null;const _=hn(f,p.x,p.y);return xe(b.data,b.width,b.height,{x:0,y:0,w:b.width,h:b.height},_,y,x,{isThreePlane:!0,threePlaneChannel:void 0,preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?oe(g==null?void 0:g.fitPoints,1,1,p.x,p.y):void 0})}}}const c=Zn(n,!1);if(t&&!t.isThreePlane&&((m=e.displaySettings)==null?void 0:m.renderMode)==="advanced-zero-dep"&&e.displaySettings.advancedZeroDep){const d=o/Math.max(1,t.width),p=a/Math.max(1,t.height);return{sourceMode:i,detectionGray:c,detectionWidth:n.width,detectionHeight:n.height,detectToMeasureX:t.width/Math.max(1,n.width),detectToMeasureY:t.height/Math.max(1,n.height),detectToDisplayX:1,detectToDisplayY:1,measureToDisplayX:d,measureToDisplayY:p,detectPointToDisplay:f=>f,measurePointToDisplay:f=>Xt(f,d,p),displayPointToDetect:f=>f,measureUsesDisplayLine:!1,measureWidth:t.width,measureHeight:t.height,refineLine:(f,y,x)=>le(c,n.width,n.height,f,y,x*.5,Math.max(4,x*.2)),measureEdge:(f,y,x,g,b)=>{const _=ao(t,f,e.displaySettings);if(!_||_.width<8||_.height<8)return null;const M=hn(y,f.x,f.y);return xe(_.data,_.width,_.height,{x:0,y:0,w:_.width,h:_.height},M,x,g,{preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?oe(b==null?void 0:b.fitPoints,t.width/Math.max(1,n.width),t.height/Math.max(1,n.height),f.x,f.y):void 0})}}}if(t&&!t.isThreePlane){const d=o/Math.max(1,t.width),p=a/Math.max(1,t.height);return{sourceMode:i,detectionGray:c,detectionWidth:n.width,detectionHeight:n.height,detectToMeasureX:t.width/Math.max(1,n.width),detectToMeasureY:t.height/Math.max(1,n.height),detectToDisplayX:1,detectToDisplayY:1,measureToDisplayX:d,measureToDisplayY:p,detectPointToDisplay:f=>f,measurePointToDisplay:f=>Xt(f,d,p),displayPointToDetect:f=>f,measureUsesDisplayLine:!1,measureWidth:t.width,measureHeight:t.height,refineLine:(f,y,x)=>le(c,n.width,n.height,f,y,x*.5,Math.max(4,x*.2)),measureEdge:(f,y,x,g,b)=>xe(t.data,t.width,t.height,f,y,x,g,{blackLevel:e.blackLevel??void 0,preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?oe(b==null?void 0:b.fitPoints,t.width/Math.max(1,n.width),t.height/Math.max(1,n.height)):void 0})}}return{sourceMode:i,detectionGray:c,detectionWidth:n.width,detectionHeight:n.height,detectToMeasureX:1,detectToMeasureY:1,detectToDisplayX:1,detectToDisplayY:1,measureToDisplayX:1,measureToDisplayY:1,detectPointToDisplay:d=>d,measurePointToDisplay:d=>d,displayPointToDetect:d=>d,measureUsesDisplayLine:!1,measureWidth:n.width,measureHeight:n.height,refineLine:(d,p,f)=>le(c,n.width,n.height,d,p,f*.5,Math.max(4,f*.2)),measureEdge:(d,p,f,y,x)=>{const g=fc(n,d);if(!g)return null;const b=hn(p,d.x,d.y);return Cl(g.data,g.width,g.height,{x:0,y:0,w:g.width,h:g.height},b,f,y,{blackLevel:e.monochromeBlackLevel??void 0,preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?oe(x==null?void 0:x.fitPoints,1,1,d.x,d.y):void 0})}}}function Ic(n,t,e){var u,h,c,m,d,p,f,y;if(!n||!t)return[];(u=e.onProgress)==null||u.call(e,"Preparing source context...",0);const i=Tc(n,t,e);if(!i)return[];(h=e.onProgress)==null||h.call(e,"Preparing source context...",.08);const r=xr(e.detectionTuning),s=Math.min(1e3,Math.max(1,e.maxRegions??1e3)),o=Math.max(4,e.maxEdges??s*4);(c=e.onProgress)==null||c.call(e,"Detecting candidates...",.12);const a=Jl(i.detectionGray,i.detectionWidth,i.detectionHeight,s,e.detectionTuning,(x,g)=>{var b;(b=e.onProgress)==null||b.call(e,x,.12+.08*Math.max(0,Math.min(1,g)))});(m=e.onProgress)==null||m.call(e,"Detecting candidates...",.2);const l=[];for(let x=0;x<a.length;x++){const g=a.length<=0?1:x/a.length;if((d=e.onProgress)==null||d.call(e,`Measuring edges: region ${x+1}/${a.length}`,.2+.72*Math.min(1,g)),l.length>=o)break;const b=a[x],_=b.corners,M=`auto-region-${x+1}`;for(let w=0;w<4&&((p=e.onProgress)==null||p.call(e,`Measuring edges: region ${x+1}/${a.length}, edge ${w+1}/4`,.2+.72*Math.min(1,(x+w/4)/Math.max(1,a.length))),!(l.length>=o));w=w+1){const P=_[w],C=_[(w+1)%4],v=C.x-P.x,k=C.y-P.y,A=Math.hypot(v,k);if(!Number.isFinite(A)||A<24)continue;const S=.125,F={x:P.x+v*S,y:P.y+k*S},T={x:C.x-v*S,y:C.y-k*S},R=Math.hypot(T.x-F.x,T.y-F.y);if(!Number.isFinite(R)||R<12)continue;const E=i.refineLine(F,T,R),O=(E!=null&&E.fitPoints?ns(E.fitPoints):null)||(E==null?void 0:E.line)||{p1:F,p2:T},U=mc(O,i.detectToMeasureX,i.detectToMeasureY),I=gc((E==null?void 0:E.fitPoints)??[],i.detectPointToDisplay),D=(I.length>=2?ns(I):null)||yc(U,i.measurePointToDisplay),B=i.measureUsesDisplayLine?D:U,V=B.p2.x-B.p1.x,Q=B.p2.y-B.p1.y,z=Math.hypot(V,Q);if(!Number.isFinite(z)||z<=1e-6)continue;const Y=D.p2.x-D.p1.x,W=D.p2.y-D.p1.y,K=Math.hypot(Y,W);if(!Number.isFinite(K)||K<=1e-6)continue;const it=!!e.distortionCurveApplied&&!!e.distortionModel,J=Y/K,ot=W/K;let Z=ot,L=-J;const G=(D.p1.x+D.p2.x)*.5,X=(D.p1.y+D.p2.y)*.5,j=i.detectPointToDisplay({x:b.centerX,y:b.centerY}),$=j.x,et=j.y;(G-$)*Z+(X-et)*L<0&&(Z=-Z,L=-L);const nt=z*.5,ct=Math.max(2,z*r.sampleHalfWidthRatio),at=Math.max(2,K*r.sampleHalfWidthRatio),rt=it?{p1:Nt(D.p1,e.distortionModel),p2:Nt(D.p2,e.distortionModel)}:void 0,Jt=rt?Math.max(1,Math.hypot(rt.p2.x-rt.p1.x,rt.p2.y-rt.p1.y)*.5):nt,Wt=it?D:U,H=it?at:ct,pt=de(Wt,H);if(!pt)continue;const mt=xt(Et(pt,2),it?n.width:i.measureWidth,it?n.height:i.measureHeight);if(!mt)continue;const gt=e.distortionCorrected&&e.distortionModel&&i.sourceMode==="rggb-raw"?yr(mt,e.distortionModel,t.width,t.height):null,ht=it?ta(e.distortionSamplingPlane??e.distortionSamplingImage??n,e.distortionModel,rt,Jt,H,!!e.sfrHasGamma,mt,e.distortionBaseImage??n,Wt):i.measureEdge(mt,B,nt,H,E);if(!ht)continue;ht.autoLikeUsed=!0;const Zt=!vl(ht);if(Zt&&!Sl(ht,e.useDeshading,0))continue;const Me=e.useNR?-1:12,Le=El([ht],Me,null,e.useDeshading,0,!0);if(!Le||Le.mtf50===null||Zt&&!Pl(Le.lsfCropped))continue;const Rt=rt?Sc(rt,e.distortionModel,Math.max(21,Math.round(K*.5))):ht.quadraticProjectionUsed?_l(I,D,Math.max(21,Math.round(K*.5))):void 0,Ee=rt&&Rt&&Rt.length>=2?Et(Rt,at+2):null,Ue=Ee?Mc(Ee):de(D,at);if(!Ue)continue;const Mn=e.distortionCorrected?mt:Ee??Et(Ue,2);let Ln={x:G+Z*(at+12),y:X+L*(at+12)},En=is(J,ot);if(Rt&&Rt.length>=3){const Se=Math.floor(Rt.length/2),De=Rt[Math.max(0,Se-1)],Be=Rt[Math.min(Rt.length-1,Se+1)],pe=Rt[Se],Sn=Be.x-De.x,Pn=Be.y-De.y,te=Math.hypot(Sn,Pn);if(te>1e-6){const Oe=Pn/te,ze=-Sn/te;En=is(Sn/te,Pn/te);const Ve={x:pe.x-$,y:pe.y-et},Ge=Ve.x*Oe+Ve.y*ze>=0?1:-1;Ln={x:pe.x+Oe*Ge*(at+12),y:pe.y+ze*Ge*(at+12)}}}l.push({id:`${M}-edge-${w+1}`,regionId:M,sourceMode:i.sourceMode,edgeIndex:w,label:Le.mtf50.toFixed(3),mtf50:Le.mtf50,angle:En,orientation:ht.orientation,edgeData:ht,sourceRect:Mn,rawSourceRect:(i.sourceMode==="rggb-raw"?gt??mt:gt)??void 0,quad:Ue,line:D,originalLine:U,curveBaseLine:rt,curvePoints:Rt,labelPoint:Ln,ridgePoints:I,outerSideMeans:b.outerSideMeans,outerSideQuads:b.outerSideQuads,distortionCorrected:e.distortionCorrected??!1})}}return(f=e.onProgress)==null||f.call(e,"Finalizing results...",.98),(y=e.onProgress)==null||y.call(e,"Finalizing results...",1),l}const Nc=n=>!n.blackLevels||n.blackLevels.length<4?null:[Number(n.blackLevels[0])||0,Number(n.blackLevels[1])||0,Number(n.blackLevels[2])||0,Number(n.blackLevels[3])||0],Wi=(n,t)=>{t instanceof ArrayBuffer&&(n.includes(t)||n.push(t))};self.onmessage=async n=>{var o,a;const{id:t,buffer:e,detect:i,options:r}=n.data,s=performance.now();try{const l=performance.now(),u=await so(e),h=performance.now()-l;let c=0,m=[];if(i&&!u.isXTrans){const p=u.isThreePlane?"three-plane":"rggb-raw",f=performance.now();m=Ic({width:u.width,height:u.height},u,{...r,sourceMode:p,forceRenderedMeasurement:!1,blackLevel:(r==null?void 0:r.blackLevel)??Nc(u),onProgress:(y,x)=>{self.postMessage({id:t,type:"progress",stage:y,progress:x})}}),c=performance.now()-f}const d=[];Wi(d,e),Wi(d,(o=u.data)==null?void 0:o.buffer),Wi(d,(a=u.floatData)==null?void 0:a.buffer),self.postMessage({id:t,type:"result",success:!0,raw:u,rawFileBuffer:e,measurements:m,timings:{decodeMs:h,detectMs:c,totalMs:performance.now()-s}},d)}catch(l){self.postMessage({id:t,type:"result",success:!1,error:(l==null?void 0:l.message)||String(l)})}};
