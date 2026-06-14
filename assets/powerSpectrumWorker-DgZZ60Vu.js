var wn=Object.defineProperty;var bn=(i,e,t)=>e in i?wn(i,e,{enumerable:!0,configurable:!0,writable:!0,value:t}):i[e]=t;var E=(i,e,t)=>bn(i,typeof e!="symbol"?e+"":e,t);class Fe{constructor(e){E(this,"size");E(this,"isPowerOfTwo");E(this,"_real");E(this,"_imag");E(this,"_scratch",null);E(this,"_rev",null);E(this,"_m",0);E(this,"_internalFFT",null);E(this,"_chirpReal",null);E(this,"_chirpImag",null);E(this,"_bReal",null);E(this,"_bImag",null);E(this,"_hanning",null);E(this,"_windowSumSq",0);this.size=e,this.isPowerOfTwo=(e&e-1)===0&&e>0,this._real=new Float32Array(e),this._imag=new Float32Array(e),this.isPowerOfTwo?this.initRadix2():this.initBluestein()}initRadix2(){const e=this.size,t=Math.log2(e);this._rev=new Uint32Array(e);for(let n=0;n<e;n++){let r=0,s=n;for(let a=0;a<t;a++)r=r<<1|s&1,s>>>=1;this._rev[n]=r}}initBluestein(){const e=this.size;this._m=Math.pow(2,Math.ceil(Math.log2(2*e-1))),this._internalFFT=new Fe(this._m),this._chirpReal=new Float32Array(e),this._chirpImag=new Float32Array(e);for(let r=0;r<e;r++){const s=-Math.PI*(r*r)/e;this._chirpReal[r]=Math.cos(s),this._chirpImag[r]=Math.sin(s)}const t=new Float32Array(this._m),n=new Float32Array(this._m);for(let r=0;r<e;r++)t[r]=this._chirpReal[r],n[r]=-this._chirpImag[r];for(let r=1;r<e;r++)t[this._m-r]=t[r],n[this._m-r]=n[r];this._internalFFT.transform(t,n),this._bReal=new Float32Array(this._internalFFT._real),this._bImag=new Float32Array(this._internalFFT._imag)}initHanning(){if(this._hanning)return;const e=this.size;this._hanning=new Float32Array(e);let t=0;for(let n=0;n<e;n++){const r=.5*(1-Math.cos(2*Math.PI*n/(e-1)));this._hanning[n]=r,t+=r*r}this._windowSumSq=t}transform(e,t){this.isPowerOfTwo?this.transformRadix2(e,t):this.transformBluestein(e,t)}transformRadix2(e,t){const n=this.size,r=this._rev,s=this._real,a=this._imag;if(e===s)for(let o=0;o<n;o++){const d=r[o];if(o<d){const u=s[o],c=a[o];s[o]=s[d],a[o]=a[d],s[d]=u,a[d]=c}}else for(let o=0;o<n;o++){const d=r[o];s[o]=e[d],a[o]=t?t[d]:0}for(let o=2;o<=n;o*=2){const d=o/2,u=-2*Math.PI/o,c=Math.cos(u),l=Math.sin(u);for(let p=0;p<n;p+=o){let h=1,m=0;for(let _=0;_<d;_++){const b=p+_,w=p+_+d,C=h*s[w]-m*a[w],y=h*a[w]+m*s[w],k=s[b],S=a[b];s[b]=k+C,a[b]=S+y,s[w]=k-C,a[w]=S-y;const A=h*c-m*l,F=h*l+m*c;h=A,m=F}}}}transformBluestein(e,t){const n=this.size,r=this._m,s=this._internalFFT,a=s._real,o=s._imag;a.fill(0),o.fill(0);for(let l=0;l<n;l++){const p=e[l],h=t?t[l]:0,m=this._chirpReal[l],_=this._chirpImag[l];a[l]=p*m-h*_,o[l]=p*_+h*m}s.transformRadix2(a,o);for(let l=0;l<r;l++){const p=s._real[l],h=s._imag[l],m=this._bReal[l],_=this._bImag[l];s._real[l]=p*m-h*_,s._imag[l]=p*_+h*m}const d=s._real,u=s._imag;for(let l=0;l<r;l++)u[l]=-u[l];s.transformRadix2(d,u);const c=1/r;for(let l=0;l<n;l++){const p=s._real[l]*c,h=-s._imag[l]*c,m=this._chirpReal[l],_=this._chirpImag[l];this._real[l]=p*m-h*_,this._imag[l]=p*_+h*m}}calculateSpectrum(e,t,n=!1){const r=this.size;let s=0;for(let c=0;c<r;c++)s+=e[c];const a=s/r;this._scratch||(this._scratch=new Float32Array(r));const o=this._scratch;if(n){this.initHanning();const c=this._hanning;for(let l=0;l<r;l++)o[l]=(e[l]-a)*c[l]}else for(let c=0;c<r;c++)o[c]=e[c]-a;this.transform(o);const d=t.length;let u=1/r;n&&this._windowSumSq>0&&(u=1/this._windowSumSq);for(let c=0;c<d;c++){const l=this._real[c],p=this._imag[c];t[c]+=(l*l+p*p)*u}}calculateSpectrumWindow(e,t,n,r=!1){const s=this.size;let a=0;for(let l=0;l<s;l++)a+=e[t+l];const o=a/s;this._scratch||(this._scratch=new Float32Array(s));const d=this._scratch;if(r){this.initHanning();const l=this._hanning;for(let p=0;p<s;p++)d[p]=(e[t+p]-o)*l[p]}else for(let l=0;l<s;l++)d[l]=e[t+l]-o;this.transform(d);const u=n.length;let c=1/s;r&&this._windowSumSq>0&&(c=1/this._windowSumSq);for(let l=0;l<u;l++){const p=this._real[l],h=this._imag[l];n[l]+=(p*p+h*h)*c}}}const xn={"Sony ILCE-7RM5":"0.82 -0.2976 -0.0719 -0.4296 1.2053 0.2532 -0.0429 0.1282 0.5774"};let Je=null;async function kn(i){return Je||(Je=(async()=>{if(typeof window.loadPyodide!="function")throw new Error("Pyodide missing: window.loadPyodide not found.");const e=await window.loadPyodide();return await e.loadPackage("numpy"),e})()),Je}var Sn=`#!/usr/bin/env python3
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
`,vn=`#!/usr/bin/env python3
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
`,Cn=`#!/usr/bin/env python3
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
`,An="/assets/joraw-BAR7UfD7.js",Pn="/assets/joraw-CZg3hVRg.wasm";const vt=512,In=xn["Sony ILCE-7RM5"].split(/\s+/).map(Number).filter(Number.isFinite);let Qe=null,Ze=null;function Tn(i){if(i.byteLength<8)return null;const e=i.getUint16(0,!1);return e===18761?!0:e===19789?!1:null}function qt(i,e,t){const n=Math.min(i.length,e+t);let r="";for(let s=e;s<n;s++){const a=i[s];if(a===0)break;r+=String.fromCharCode(a)}return r.trim()}function H(i,e,t,n,r){const s=t===1||t===2||t===7?1:t===3||t===8?2:t===4||t===9?4:0;if(!s)return[];const a=s*n,o=a<=4?r:i.getUint32(r,e);if(o<0||o+a>i.byteLength)return[];const d=[];for(let u=0;u<n;u++){const c=o+u*s;t===1||t===2||t===7?d.push(i.getUint8(c)):t===3?d.push(i.getUint16(c,e)):t===8?d.push(i.getInt16(c,e)):t===4?d.push(i.getUint32(c,e)):t===9&&d.push(i.getInt32(c,e))}return d}function Ct(i,e,t,n,r,s){if(n!==2||r<=0)return"";const a=r<=4?s:e.getUint32(s,t);return a<0||a>=i.length?"":qt(i,a,r)}function Yt(i){const e=new Uint8Array(i),t=new DataView(i),n=Tn(t);if(n===null||t.getUint16(2,n)!==42)return null;const s=l=>t.getUint16(l,n),a=l=>t.getUint32(l,n),o=[a(4)],d=new Set;let u="",c="";for(;o.length;){const l=o.pop();if(d.has(l)||l<=0||l+2>t.byteLength)continue;d.add(l);const p=s(l);if(l+2+p*12+4>t.byteLength)continue;const h=new Map;for(let k=0;k<p;k++){const S=l+2+k*12,A=s(S),F=s(S+2),X=a(S+4),x=S+8;h.set(A,{type:F,count:X,valueOffset:x})}const m=h.get(271),_=h.get(272);m&&!u&&(u=Ct(e,t,n,m.type,m.count,m.valueOffset)),_&&!c&&(c=Ct(e,t,n,_.type,_.count,_.valueOffset));const b=h.get(330);if(b){const k=H(t,n,b.type,b.count,b.valueOffset);for(const S of k)o.push(S)}const w=h.get(259),C=h.get(262);if(w&&C){const k=H(t,n,w.type,w.count,w.valueOffset)[0],S=H(t,n,C.type,C.count,C.valueOffset)[0];if(k===32766&&S===32803){const A=H(t,n,h.get(256).type,h.get(256).count,h.get(256).valueOffset)[0],F=H(t,n,h.get(257).type,h.get(257).count,h.get(257).valueOffset)[0],X=H(t,n,h.get(258).type,h.get(258).count,h.get(258).valueOffset)[0],x=H(t,n,h.get(273).type,h.get(273).count,h.get(273).valueOffset)[0],T=H(t,n,h.get(279).type,h.get(279).count,h.get(279).valueOffset)[0],V=h.get(33422)?H(t,n,h.get(33422).type,h.get(33422).count,h.get(33422).valueOffset):[0,1,1,2];h.get(29456)&&H(t,n,h.get(29456).type,h.get(29456).count,h.get(29456).valueOffset);const z=h.get(50717)?H(t,n,h.get(50717).type,h.get(50717).count,h.get(50717).valueOffset)[0]:16383,K=h.get(50719)?H(t,n,h.get(50719).type,h.get(50719).count,h.get(50719).valueOffset):[],G=h.get(50720)?H(t,n,h.get(50720).type,h.get(50720).count,h.get(50720).valueOffset):[];if(x+vt+16>e.length||x+T>e.length)return null;const f=x+vt,L=qt(e,f,4),U=e[f+8]<<8|e[f+9],M=e[f+10]<<8|e[f+11],P=e[f+12]<<8|e[f+13],q=e[f+14]<<8|e[f+15],qe=P>>4&63,pn=q>>13,fn=q>>10&3,Ie=M*2,se=T>=4?(e[x]|e[x+1]<<8|e[x+2]<<16|e[x+3]<<24)>>>0:0,gn=U===A&&Ie===F;let kt=!1;if(se>=1&&se<=16&&T>=8+se*24){const ae=new Map,Ye=new Map;let ke=!0;for(let Se=0;Se<se;Se++){const Q=x+8+Se*24,oe=t.getUint32(Q+8,!0),le=t.getUint32(Q+12,!0),$e=t.getUint32(Q+16,!0),Te=t.getUint32(Q+20,!0);if(!$e||!Te||oe+$e>A||le+Te>F){ke=!1;break}const St=ae.get(le);if(St!==void 0&&St!==Te){ke=!1;break}ae.set(le,Te),Ye.set(le,(Ye.get(le)||0)+$e)}if(ke){const Se=Array.from(ae.keys()).sort((oe,le)=>oe-le);let Q=0;for(const oe of Se){if(oe!==Q||Ye.get(oe)!==A){ke=!1;break}Q+=ae.get(oe)}kt=ke&&Q===F}}const mn=se>=1&&se<=16&&U>0&&Ie>0&&A%U===0&&F%Ie===0&&se===A/U*(F/Ie);if(L!=="A000"&&L!=="0000"||!gn&&!mn&&!kt||qe!==16||pn!==3||fn!==3)return null;const _n=[1024,1024,1024,1024],yn=V.slice(0,4).map(ae=>ae===0?"R":ae===2?"B":"G").join("")||"RGGB";return{width:A,height:F,bitsPerSample:X,compression:k,photometric:S,blackLevel:_n,whiteLevel:Number(z||16383),cfaPattern:yn,defaultCropOrigin:K.length>=2?[Number(K[0]),Number(K[1])]:void 0,defaultCropSize:G.length>=2?[Number(G[0]),Number(G[1])]:void 0,make:u||"SONY",model:c||"ILCE-7M5"}}}const y=a(l+2+p*12);y&&o.push(y)}return null}async function Ln(i){return Qe||(Qe=(async()=>{const e=await kn();return e.__jtrSonyCrawHqDecoderReady||(await e.FS.mkdirTree("/sony_craw_hq"),await e.FS.writeFile("/sony_craw_hq/llvc3_bitstream_probe.py",Sn),await e.FS.writeFile("/sony_craw_hq/llvc3_entropy.py",vn),await e.FS.writeFile("/sony_craw_hq/llvc3_math.py",Cn),await e.runPythonAsync(`
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
`),e.__jtrSonyCrawHqDecoderReady=!0),e})()),Qe}function Rn(i){return Yt(i)}async function On(){return Ze||(Ze=import(An).then(i=>{const e=i.default;if(typeof e!="function")throw new Error("Sony cRAW HQ LibRaw WASM module is missing its initializer");return e({locateFile:(t,n)=>t.endsWith("joraw.wasm")?Pn:n+t})})),Ze}function Mn(i){return i==="ILCE-7M5"?In:null}function $t(i,e,t,n){var d;if(i.length!==e.width*e.height)throw new Error(`Sony cRAW HQ decoded size mismatch: got ${i.length}, expected ${e.width*e.height}`);const r=e.model||"ILCE-7M5",s=r.startsWith("Sony ")?r:`Sony ${r}`,a=Mn(r),o={...n||{},make:e.make||(n==null?void 0:n.camera_make)||"SONY",model:r,camera_make:e.make||(n==null?void 0:n.camera_make)||"SONY",camera_model:r,UniqueCameraModel:s,sourceFormat:t==="libraw-wasm"?"Sony cRAW HQ / LLVC3 (LibRaw WASM)":"Sony cRAW HQ / LLVC3",sonyCrawHq:{...e,decodeBackend:t},color_desc:e.cfaPattern,black_level_per_channel:e.blackLevel,white_level:e.whiteLevel,color_matrix:a&&a.length===9?a:void 0,idata:{filters:2492765332,colors:3},color_data:{...(n==null?void 0:n.color_data)||{},black:1024,cblack_rawpy_style:e.blackLevel,dng_levels:{...((d=n==null?void 0:n.color_data)==null?void 0:d.dng_levels)||{},dng_cblack:e.blackLevel,dng_whitelevel:e.whiteLevel}}};return{data:i,width:e.width,height:e.height,bayerPattern:e.cfaPattern,blackLevels:e.blackLevel,whiteLevel:e.whiteLevel,metadata:o,isThreePlane:!1,isXTrans:!1}}async function Dn(i,e,t){const n=typeof performance<"u"?performance.now():Date.now(),r=await On(),s=typeof performance<"u"?performance.now():Date.now(),a=r.LibRaw||r.JoRaw;if(!a)throw new Error("Sony cRAW HQ LibRaw WASM class not found");const o=new a;try{const d=new Uint8Array(i);o.open(d,{});const u=typeof performance<"u"?performance.now():Date.now();let c=null;try{c=o.metadata(!0)}catch(w){console.warn("[Sony cRAW HQ] fast WASM metadata read failed",w)}const l=typeof performance<"u"?performance.now():Date.now(),p=o.getRawImage(),h=typeof performance<"u"?performance.now():Date.now();if(!p||!p.data)throw new Error("Sony cRAW HQ LibRaw WASM returned no raw image");const m=p.data instanceof Uint16Array?p.data:new Uint16Array(p.data.buffer,p.data.byteOffset||0,p.data.byteLength/2),_=typeof performance<"u"?performance.now():Date.now();if(p.width!==e.width||p.height!==e.height)throw new Error(`Sony cRAW HQ LibRaw WASM dimensions mismatch: got ${p.width}x${p.height}, expected ${e.width}x${e.height}`);const b=typeof performance<"u"?performance.now():Date.now();return console.info("[Sony cRAW HQ] fast decode timings",{width:e.width,height:e.height,backend:"libraw-wasm",wasmReadyMs:Math.round(s-n),openMs:Math.round(u-s),metadataMs:Math.round(l-u),unpackMs:Math.round(h-l),copyMs:Math.round(_-h),totalMs:Math.round(b-n)}),{rawImageData:$t(m,e,"libraw-wasm",c),info:e}}finally{typeof o.delete=="function"?o.delete():typeof o.close=="function"&&o.close()}}async function Fn(i,e,t){const n=typeof performance<"u"?performance.now():Date.now(),r=await Ln(),s=typeof performance<"u"?performance.now():Date.now(),a=new Uint8Array(i),o=await fetch(new URL("/assets/sony_llvc3_static_lut4096_padded_u16-FsVBk-IV.bin",import.meta.url));if(!o.ok)throw new Error(`Failed to load Sony LLVC3 sample LUT: HTTP ${o.status}`);const d=new Uint8Array(await o.arrayBuffer()),u=typeof performance<"u"?performance.now():Date.now();r.globals.set("jtr_sony_arw_bytes",a),r.globals.set("jtr_sony_lut_bytes",d);const c=await r.runPythonAsync("jtr_decode_sony_craw_hq(jtr_sony_arw_bytes.to_py(), jtr_sony_lut_bytes.to_py())"),l=typeof performance<"u"?performance.now():Date.now(),p=c.toJs();typeof c.destroy=="function"&&c.destroy(),r.globals.delete("jtr_sony_arw_bytes"),r.globals.delete("jtr_sony_lut_bytes");const h=new Uint8Array(p.byteLength);h.set(p);const m=new Uint16Array(h.buffer),_=typeof performance<"u"?performance.now():Date.now(),b=typeof performance<"u"?performance.now():Date.now();return console.info("[Sony cRAW HQ] decode timings",{width:e.width,height:e.height,backend:"pyodide",pyodideReadyMs:Math.round(s-n),lutLoadMs:Math.round(u-s),llvc3DecodeMs:Math.round(l-u),copyMs:Math.round(_-l),totalMs:Math.round(b-n)}),{rawImageData:$t(m,e,"pyodide"),info:e}}async function En(i,e){const t=Yt(i);if(!t)return null;try{return await Dn(i,t,e)}catch(n){return console.warn("[Sony cRAW HQ] fast WASM decode failed; falling back to Pyodide",n),Fn(i,t)}}async function Nn(i,e){return En(i,e)}const Un={1:1,2:1,3:2,4:4,5:8,7:1,9:4,10:8,11:4,12:8},Bn=34713,lt=155;let et=null;function At(i,e,t){return ye(i,e,t,Math.min(t.count,256)).map(r=>r&255).filter(r=>r!==0).map(r=>String.fromCharCode(r)).join("").trim()}function ye(i,e,t,n=Number.MAX_SAFE_INTEGER){const r=Un[t.type]||0;if(!r)return[];const a=t.count*r<=4?t.valueOffset:i.getUint32(t.valueOffset,e),o=Math.min(t.count,n),d=[];for(let u=0;u<o;u++){const c=a+u*r;if(c<0||c+r>i.byteLength)break;switch(t.type){case 1:case 2:case 7:d.push(i.getUint8(c));break;case 3:d.push(i.getUint16(c,e));break;case 4:d.push(i.getUint32(c,e));break;case 9:d.push(i.getInt32(c,e));break;case 5:{const l=i.getUint32(c,e),p=i.getUint32(c+4,e);d.push(p?l/p:0);break}case 10:{const l=i.getInt32(c,e),p=i.getInt32(c+4,e);d.push(p?l/p:0);break}}}return d}function Pt(i,e,t){if(t<=0||t+2>i.byteLength)return null;const n=i.getUint16(t,e),r=2+n*12+4;if(n>4096||t+r>i.byteLength)return null;const s=new Map;for(let a=0;a<n;a++){const o=t+2+a*12,d=i.getUint16(o,e);s.set(d,{tag:d,type:i.getUint16(o+2,e),count:i.getUint32(o+4,e),valueOffset:o+8})}return s}function Vn(i,e,t){if(t<=0||t+2>i.byteLength)return 0;const n=i.getUint16(t,e),r=t+2+n*12;return r+4>i.byteLength?0:i.getUint32(r,e)}function zn(i){const e=["R","G","B","C","M","Y","W"];return i.length<4?"RGGB":i.slice(0,4).map(t=>e[t]||"G").join("")}function de(i,e,t,n){const r=i.get(n);if(!r)return null;const s=ye(e,t,r,1);return s.length?s[0]:null}function Jt(i){if(i.byteLength<16)return null;const e=new DataView(i),t=e.getUint16(0,!1),n=t===18761;if(!n&&t!==19789||e.getUint16(2,n)!==42)return null;const r=e.getUint32(4,n),s=Pt(e,n,r);if(!s)return null;const a=s.get(271),o=s.get(272),d=a?At(e,n,a):"",u=o?At(e,n,o):"";if(!/NIKON/i.test(d+" "+u))return null;const c=de(s,e,n,274)||1,l=s.get(330);if(!l)return null;const p=ye(e,n,l,16);for(const m of p){const _=Pt(e,n,m);if(!_)continue;const b=de(_,e,n,259)||0,w=de(_,e,n,273)||0,C=de(_,e,n,279)||0,y=de(_,e,n,256)||0,k=de(_,e,n,257)||0;if(b!==Bn||!w||!C||y<=0||k<=0||w+C>i.byteLength)continue;const S=de(_,e,n,258)||14,A=_.get(33422),F=A?zn(ye(e,n,A,4)):"RGGB",X=_.get(50714),x=_.get(50717),T=X?ye(e,n,X,4):[],V=x?ye(e,n,x,1):[],z=T.length&&Number(T[0])||0,K=[Number(T[0])||z,Number(T[1])||z,Number(T[2])||z,Number(T[3])||z],G=w+lt,f=G+3<i.byteLength?e.getUint8(G+3):0;return{make:d,model:u,width:y,height:k,stripOffset:w,stripByteCount:C,bayerPattern:F,blackLevels:K,whiteLevel:Math.max(1,Number(V[0])||(1<<Math.min(S,15))-1||16383),orientation:c,compression:b,firstBp:f,isHeStar:f>=1&&f<=2}}return Vn(e,n,r)>0,null}function Hn(){return{wasi_snapshot_preview1:{fd_close:()=>0,fd_seek:(i,e,t,n,r)=>0,fd_write:()=>0,proc_exit:i=>{throw new Error(`Nikon HE WASM exited with code ${i}`)}}}}async function Gn(){return et||(et=(async()=>{const i=new URL("/assets/nikon-he-decoder-DUHNoGN2.wasm",import.meta.url),e=await fetch(i);if(!e.ok)throw new Error(`Failed to load Nikon HE WASM: ${e.status}`);const t=await e.arrayBuffer(),r=(await WebAssembly.instantiate(t,Hn())).instance.exports;if(!r.memory||!r.nikon_he_decode||!r.nikon_he_malloc||!r.nikon_he_free)throw new Error("Nikon HE WASM exports are incomplete.");return r})()),et}function Wn(i){return Jt(i)}async function jn(i){const e=Jt(i);if(!e)return null;const t=e.stripOffset+lt,n=e.stripByteCount-lt;if(n<=12||t+n>i.byteLength)throw new Error("Invalid Nikon HE precinct bounds.");const r=await Gn(),s=e.width*e.height*2,a=r.nikon_he_malloc(n),o=r.nikon_he_malloc(s);try{new Uint8Array(r.memory.buffer).set(new Uint8Array(i,t,n),a);const u=r.nikon_he_decode(a,n,e.width,e.height,o);if(u<0)throw new Error(`Nikon HE decoder failed with code ${u}.`);const c=new Uint16Array(r.memory.buffer,o,e.width*e.height),l=new Uint16Array(c);let p=65535,h=0;const m=Math.max(1,Math.floor(l.length/2e5));for(let _=0;_<l.length;_+=m){const b=l[_];b<p&&(p=b),b>h&&(h=b)}if(p===h)throw new Error("Nikon HE decoder returned a constant image; this HE variant is not supported by the current WASM core.");return{data:l,width:e.width,height:e.height,bayerPattern:e.bayerPattern,blackLevels:e.blackLevels,whiteLevel:e.whiteLevel,metadata:{description:`Nikon ${e.isHeStar?"HE*":"HE"} RAW decoded by standalone WASM`,camera_make:e.make,camera_model:e.model,make:e.make,model:e.model,orientation:e.orientation,compression:e.compression,raw_width:e.width,raw_height:e.height,width:e.width,height:e.height,color_desc:e.bayerPattern,white_level:e.whiteLevel,black_level_per_channel:e.blackLevels,nikon_he:{isHeStar:e.isHeStar,firstBp:e.firstBp,stripOffset:e.stripOffset,stripByteCount:e.stripByteCount,precinctOffset:t,precinctSize:n}},isThreePlane:!1}}finally{r.nikon_he_free(o),r.nikon_he_free(a)}}var He=typeof self<"u"?self:global;const Pe=typeof navigator<"u",Xn=Pe&&typeof HTMLImageElement>"u",Ee=!(typeof global>"u"||typeof process>"u"||!process.versions||!process.versions.node),Ge=He.Buffer,Le=He.BigInt,We=!!Ge,Kn=i=>i;function Ne(i,e=Kn){if(Ee)try{return typeof require=="function"?Promise.resolve(e(require(i))):import(i).then(e)}catch{console.warn(`Couldn't load ${i}`)}}let ft=He.fetch;const qn=i=>ft=i;if(!He.fetch){const i=Ne("http",n=>n),e=Ne("https",n=>n),t=(n,{headers:r}={})=>new Promise(async(s,a)=>{let{port:o,hostname:d,pathname:u,protocol:c,search:l}=new URL(n);const p={method:"GET",hostname:d,path:encodeURI(u)+l,headers:r};o!==""&&(p.port=Number(o));const h=(c==="https:"?await e:await i).request(p,m=>{if(m.statusCode===301||m.statusCode===302){let _=new URL(m.headers.location,n).toString();return t(_,{headers:r}).then(s).catch(a)}s({status:m.statusCode,arrayBuffer:()=>new Promise(_=>{let b=[];m.on("data",w=>b.push(w)),m.on("end",()=>_(Buffer.concat(b)))})})});h.on("error",a),h.end()});qn(t)}function g(i,e,t){return e in i?Object.defineProperty(i,e,{value:t,enumerable:!0,configurable:!0,writable:!0}):i[e]=t,i}const Ue=i=>Qt(i)?void 0:i,Yn=i=>i!==void 0;function Qt(i){return i===void 0||(i instanceof Map?i.size===0:Object.values(i).filter(Yn).length===0)}function O(i){let e=new Error(i);throw delete e.stack,e}function ce(i){return(i=function(e){for(;e.endsWith("\0");)e=e.slice(0,-1);return e}(i).trim())===""?void 0:i}function dt(i){let e=function(t){let n=0;return t.ifd0.enabled&&(n+=1024),t.exif.enabled&&(n+=2048),t.makerNote&&(n+=2048),t.userComment&&(n+=1024),t.gps.enabled&&(n+=512),t.interop.enabled&&(n+=100),t.ifd1.enabled&&(n+=1024),n+2048}(i);return i.jfif.enabled&&(e+=50),i.xmp.enabled&&(e+=2e4),i.iptc.enabled&&(e+=14e3),i.icc.enabled&&(e+=6e3),e}const ct=i=>String.fromCharCode.apply(null,i),It=typeof TextDecoder<"u"?new TextDecoder("utf-8"):void 0;function Zt(i){return It?It.decode(i):We?Buffer.from(i).toString("utf8"):decodeURIComponent(escape(ct(i)))}class B{static from(e,t){return e instanceof this&&e.le===t?e:new B(e,void 0,void 0,t)}constructor(e,t=0,n,r){if(typeof r=="boolean"&&(this.le=r),Array.isArray(e)&&(e=new Uint8Array(e)),e===0)this.byteOffset=0,this.byteLength=0;else if(e instanceof ArrayBuffer){n===void 0&&(n=e.byteLength-t);let s=new DataView(e,t,n);this._swapDataView(s)}else if(e instanceof Uint8Array||e instanceof DataView||e instanceof B){n===void 0&&(n=e.byteLength-t),(t+=e.byteOffset)+n>e.byteOffset+e.byteLength&&O("Creating view outside of available memory in ArrayBuffer");let s=new DataView(e.buffer,t,n);this._swapDataView(s)}else if(typeof e=="number"){let s=new DataView(new ArrayBuffer(e));this._swapDataView(s)}else O("Invalid input argument for BufferView: "+e)}_swapArrayBuffer(e){this._swapDataView(new DataView(e))}_swapBuffer(e){this._swapDataView(new DataView(e.buffer,e.byteOffset,e.byteLength))}_swapDataView(e){this.dataView=e,this.buffer=e.buffer,this.byteOffset=e.byteOffset,this.byteLength=e.byteLength}_lengthToEnd(e){return this.byteLength-e}set(e,t,n=B){return e instanceof DataView||e instanceof B?e=new Uint8Array(e.buffer,e.byteOffset,e.byteLength):e instanceof ArrayBuffer&&(e=new Uint8Array(e)),e instanceof Uint8Array||O("BufferView.set(): Invalid data argument."),this.toUint8().set(e,t),new n(this,t,e.byteLength)}subarray(e,t){return t=t||this._lengthToEnd(e),new B(this,e,t)}toUint8(){return new Uint8Array(this.buffer,this.byteOffset,this.byteLength)}getUint8Array(e,t){return new Uint8Array(this.buffer,this.byteOffset+e,t)}getString(e=0,t=this.byteLength){return Zt(this.getUint8Array(e,t))}getLatin1String(e=0,t=this.byteLength){let n=this.getUint8Array(e,t);return ct(n)}getUnicodeString(e=0,t=this.byteLength){const n=[];for(let r=0;r<t&&e+r<this.byteLength;r+=2)n.push(this.getUint16(e+r));return ct(n)}getInt8(e){return this.dataView.getInt8(e)}getUint8(e){return this.dataView.getUint8(e)}getInt16(e,t=this.le){return this.dataView.getInt16(e,t)}getInt32(e,t=this.le){return this.dataView.getInt32(e,t)}getUint16(e,t=this.le){return this.dataView.getUint16(e,t)}getUint32(e,t=this.le){return this.dataView.getUint32(e,t)}getFloat32(e,t=this.le){return this.dataView.getFloat32(e,t)}getFloat64(e,t=this.le){return this.dataView.getFloat64(e,t)}getFloat(e,t=this.le){return this.dataView.getFloat32(e,t)}getDouble(e,t=this.le){return this.dataView.getFloat64(e,t)}getUintBytes(e,t,n){switch(t){case 1:return this.getUint8(e,n);case 2:return this.getUint16(e,n);case 4:return this.getUint32(e,n);case 8:return this.getUint64&&this.getUint64(e,n)}}getUint(e,t,n){switch(t){case 8:return this.getUint8(e,n);case 16:return this.getUint16(e,n);case 32:return this.getUint32(e,n);case 64:return this.getUint64&&this.getUint64(e,n)}}toString(e){return this.dataView.toString(e,this.constructor.name)}ensureChunk(){}}function ht(i,e){O(`${i} '${e}' was not loaded, try using full build of exifr.`)}class gt extends Map{constructor(e){super(),this.kind=e}get(e,t){return this.has(e)||ht(this.kind,e),t&&(e in t||function(n,r){O(`Unknown ${n} '${r}'.`)}(this.kind,e),t[e].enabled||ht(this.kind,e)),super.get(e)}keyList(){return Array.from(this.keys())}}var Y=new gt("file parser"),R=new gt("segment parser"),$=new gt("file reader");function $n(i,e){return typeof i=="string"?Tt(i,e):Pe&&!Xn&&i instanceof HTMLImageElement?Tt(i.src,e):i instanceof Uint8Array||i instanceof ArrayBuffer||i instanceof DataView?new B(i):Pe&&i instanceof Blob?ut(i,e,"blob",fe):void O("Invalid input argument")}function Tt(i,e){return(t=i).startsWith("data:")||t.length>1e4?pt(i,e,"base64"):Ee&&i.includes("://")?ut(i,e,"url",pe):Ee?pt(i,e,"fs"):Pe?ut(i,e,"url",pe):void O("Invalid input argument");var t}async function ut(i,e,t,n){return $.has(t)?pt(i,e,t):n?async function(r,s){let a=await s(r);return new B(a)}(i,n):void O(`Parser ${t} is not loaded`)}async function pt(i,e,t){let n=new($.get(t))(i,e);return await n.read(),n}const pe=i=>ft(i).then(e=>e.arrayBuffer()),fe=i=>new Promise((e,t)=>{let n=new FileReader;n.onloadend=()=>e(n.result||new ArrayBuffer),n.onerror=t,n.readAsArrayBuffer(i)});class Jn extends Map{get tagKeys(){return this.allKeys||(this.allKeys=Array.from(this.keys())),this.allKeys}get tagValues(){return this.allValues||(this.allValues=Array.from(this.values())),this.allValues}}function I(i,e,t){let n=new Jn;for(let[r,s]of t)n.set(r,s);if(Array.isArray(e))for(let r of e)i.set(r,n);else i.set(e,n);return n}function ge(i,e,t){let n,r=i.get(e);for(n of t)r.set(n[0],n[1])}const D=new Map,j=new Map,ie=new Map,Z=["chunked","firstChunkSize","firstChunkSizeNode","firstChunkSizeBrowser","chunkSize","chunkLimit"],be=["jfif","xmp","icc","iptc","ihdr"],me=["tiff",...be],v=["ifd0","ifd1","exif","gps","interop"],ee=[...me,...v],te=["makerNote","userComment"],xe=["translateKeys","translateValues","reviveValues","multiSegment"],ne=[...xe,"sanitize","mergeOutput","silentErrors"];class en{get translate(){return this.translateKeys||this.translateValues||this.reviveValues}}class ve extends en{get needed(){return this.enabled||this.deps.size>0}constructor(e,t,n,r){if(super(),g(this,"enabled",!1),g(this,"skip",new Set),g(this,"pick",new Set),g(this,"deps",new Set),g(this,"translateKeys",!1),g(this,"translateValues",!1),g(this,"reviveValues",!1),this.key=e,this.enabled=t,this.parse=this.enabled,this.applyInheritables(r),this.canBeFiltered=v.includes(e),this.canBeFiltered&&(this.dict=D.get(e)),n!==void 0)if(Array.isArray(n))this.parse=this.enabled=!0,this.canBeFiltered&&n.length>0&&this.translateTagSet(n,this.pick);else if(typeof n=="object"){if(this.enabled=!0,this.parse=n.parse!==!1,this.canBeFiltered){let{pick:s,skip:a}=n;s&&s.length>0&&this.translateTagSet(s,this.pick),a&&a.length>0&&this.translateTagSet(a,this.skip)}this.applyInheritables(n)}else n===!0||n===!1?this.parse=this.enabled=n:O(`Invalid options argument: ${n}`)}applyInheritables(e){let t,n;for(t of xe)n=e[t],n!==void 0&&(this[t]=n)}translateTagSet(e,t){if(this.dict){let n,r,{tagKeys:s,tagValues:a}=this.dict;for(n of e)typeof n=="string"?(r=a.indexOf(n),r===-1&&(r=s.indexOf(Number(n))),r!==-1&&t.add(Number(s[r]))):t.add(n)}else for(let n of e)t.add(n)}finalizeFilters(){!this.enabled&&this.deps.size>0?(this.enabled=!0,Be(this.pick,this.deps)):this.enabled&&this.pick.size>0&&Be(this.pick,this.deps)}}var N={jfif:!1,tiff:!0,xmp:!1,icc:!1,iptc:!1,ifd0:!0,ifd1:!1,exif:!0,gps:!0,interop:!1,ihdr:void 0,makerNote:!1,userComment:!1,multiSegment:!1,skip:[],pick:[],translateKeys:!0,translateValues:!0,reviveValues:!0,sanitize:!0,mergeOutput:!0,silentErrors:!0,chunked:!0,firstChunkSize:void 0,firstChunkSizeNode:512,firstChunkSizeBrowser:65536,chunkSize:65536,chunkLimit:5},Lt=new Map;class _e extends en{static useCached(e){let t=Lt.get(e);return t!==void 0||(t=new this(e),Lt.set(e,t)),t}constructor(e){super(),e===!0?this.setupFromTrue():e===void 0?this.setupFromUndefined():Array.isArray(e)?this.setupFromArray(e):typeof e=="object"?this.setupFromObject(e):O(`Invalid options argument ${e}`),this.firstChunkSize===void 0&&(this.firstChunkSize=Pe?this.firstChunkSizeBrowser:this.firstChunkSizeNode),this.mergeOutput&&(this.ifd1.enabled=!1),this.filterNestedSegmentTags(),this.traverseTiffDependencyTree(),this.checkLoadedPlugins()}setupFromUndefined(){let e;for(e of Z)this[e]=N[e];for(e of ne)this[e]=N[e];for(e of te)this[e]=N[e];for(e of ee)this[e]=new ve(e,N[e],void 0,this)}setupFromTrue(){let e;for(e of Z)this[e]=N[e];for(e of ne)this[e]=N[e];for(e of te)this[e]=!0;for(e of ee)this[e]=new ve(e,!0,void 0,this)}setupFromArray(e){let t;for(t of Z)this[t]=N[t];for(t of ne)this[t]=N[t];for(t of te)this[t]=N[t];for(t of ee)this[t]=new ve(t,!1,void 0,this);this.setupGlobalFilters(e,void 0,v)}setupFromObject(e){let t;for(t of(v.ifd0=v.ifd0||v.image,v.ifd1=v.ifd1||v.thumbnail,Object.assign(this,e),Z))this[t]=tt(e[t],N[t]);for(t of ne)this[t]=tt(e[t],N[t]);for(t of te)this[t]=tt(e[t],N[t]);for(t of me)this[t]=new ve(t,N[t],e[t],this);for(t of v)this[t]=new ve(t,N[t],e[t],this.tiff);this.setupGlobalFilters(e.pick,e.skip,v,ee),e.tiff===!0?this.batchEnableWithBool(v,!0):e.tiff===!1?this.batchEnableWithUserValue(v,e):Array.isArray(e.tiff)?this.setupGlobalFilters(e.tiff,void 0,v):typeof e.tiff=="object"&&this.setupGlobalFilters(e.tiff.pick,e.tiff.skip,v)}batchEnableWithBool(e,t){for(let n of e)this[n].enabled=t}batchEnableWithUserValue(e,t){for(let n of e){let r=t[n];this[n].enabled=r!==!1&&r!==void 0}}setupGlobalFilters(e,t,n,r=n){if(e&&e.length){for(let a of r)this[a].enabled=!1;let s=Rt(e,n);for(let[a,o]of s)Be(this[a].pick,o),this[a].enabled=!0}else if(t&&t.length){let s=Rt(t,n);for(let[a,o]of s)Be(this[a].skip,o)}}filterNestedSegmentTags(){let{ifd0:e,exif:t,xmp:n,iptc:r,icc:s}=this;this.makerNote?t.deps.add(37500):t.skip.add(37500),this.userComment?t.deps.add(37510):t.skip.add(37510),n.enabled||e.skip.add(700),r.enabled||e.skip.add(33723),s.enabled||e.skip.add(34675)}traverseTiffDependencyTree(){let{ifd0:e,exif:t,gps:n,interop:r}=this;r.needed&&(t.deps.add(40965),e.deps.add(40965)),t.needed&&e.deps.add(34665),n.needed&&e.deps.add(34853),this.tiff.enabled=v.some(s=>this[s].enabled===!0)||this.makerNote||this.userComment;for(let s of v)this[s].finalizeFilters()}get onlyTiff(){return!be.map(e=>this[e].enabled).some(e=>e===!0)&&this.tiff.enabled}checkLoadedPlugins(){for(let e of me)this[e].enabled&&!R.has(e)&&ht("segment parser",e)}}function Rt(i,e){let t,n,r,s,a=[];for(r of e){for(s of(t=D.get(r),n=[],t))(i.includes(s[0])||i.includes(s[1]))&&n.push(s[0]);n.length&&a.push([r,n])}return a}function tt(i,e){return i!==void 0?i:e!==void 0?e:void 0}function Be(i,e){for(let t of e)i.add(t)}g(_e,"default",N);class re{constructor(e){g(this,"parsers",{}),g(this,"output",{}),g(this,"errors",[]),g(this,"pushToErrors",t=>this.errors.push(t)),this.options=_e.useCached(e)}async read(e){this.file=await $n(e,this.options)}setup(){if(this.fileParser)return;let{file:e}=this,t=e.getUint16(0);for(let[n,r]of Y)if(r.canHandle(e,t))return this.fileParser=new r(this.options,this.file,this.parsers),e[n]=!0;this.file.close&&this.file.close(),O("Unknown file format")}async parse(){let{output:e,errors:t}=this;return this.setup(),this.options.silentErrors?(await this.executeParsers().catch(this.pushToErrors),t.push(...this.fileParser.errors)):await this.executeParsers(),this.file.close&&this.file.close(),this.options.silentErrors&&t.length>0&&(e.errors=t),Ue(e)}async executeParsers(){let{output:e}=this;await this.fileParser.parse();let t=Object.values(this.parsers).map(async n=>{let r=await n.parse();n.assignToOutput(e,r)});this.options.silentErrors&&(t=t.map(n=>n.catch(this.pushToErrors))),await Promise.all(t)}async extractThumbnail(){this.setup();let{options:e,file:t}=this,n=R.get("tiff",e);var r;if(t.tiff?r={start:0,type:"tiff"}:t.jpeg&&(r=await this.fileParser.getOrFindSegment("tiff")),r===void 0)return;let s=await this.fileParser.ensureSegmentChunk(r),a=this.parsers.tiff=new n(s,e,t),o=await a.extractThumbnail();return t.close&&t.close(),o}}async function je(i,e){let t=new re(e);return await t.read(i),t.parse()}var Qn=Object.freeze({__proto__:null,parse:je,Exifr:re,fileParsers:Y,segmentParsers:R,fileReaders:$,tagKeys:D,tagValues:j,tagRevivers:ie,createDictionary:I,extendDictionary:ge,fetchUrlAsArrayBuffer:pe,readBlobAsArrayBuffer:fe,chunkedProps:Z,otherSegments:be,segments:me,tiffBlocks:v,segmentsAndBlocks:ee,tiffExtractables:te,inheritables:xe,allFormatters:ne,Options:_e});class Xe{constructor(e,t,n){g(this,"errors",[]),g(this,"ensureSegmentChunk",async r=>{let s=r.start,a=r.size||65536;if(this.file.chunked)if(this.file.available(s,a))r.chunk=this.file.subarray(s,a);else try{r.chunk=await this.file.readChunk(s,a)}catch(o){O(`Couldn't read segment: ${JSON.stringify(r)}. ${o.message}`)}else this.file.byteLength>s+a?r.chunk=this.file.subarray(s,a):r.size===void 0?r.chunk=this.file.subarray(s):O("Segment unreachable: "+JSON.stringify(r));return r.chunk}),this.extendOptions&&this.extendOptions(e),this.options=e,this.file=t,this.parsers=n}injectSegment(e,t){this.options[e].enabled&&this.createParser(e,t)}createParser(e,t){let n=new(R.get(e))(t,this.options,this.file);return this.parsers[e]=n}createParsers(e){for(let t of e){let{type:n,chunk:r}=t,s=this.options[n];if(s&&s.enabled){let a=this.parsers[n];a&&a.append||a||this.createParser(n,r)}}}async readSegments(e){let t=e.map(this.ensureSegmentChunk);await Promise.all(t)}}class W{static findPosition(e,t){let n=e.getUint16(t+2)+2,r=typeof this.headerLength=="function"?this.headerLength(e,t,n):this.headerLength,s=t+r,a=n-r;return{offset:t,length:n,headerLength:r,start:s,size:a,end:s+a}}static parse(e,t={}){return new this(e,new _e({[this.type]:t}),e).parse()}normalizeInput(e){return e instanceof B?e:new B(e)}constructor(e,t={},n){g(this,"errors",[]),g(this,"raw",new Map),g(this,"handleError",r=>{if(!this.options.silentErrors)throw r;this.errors.push(r.message)}),this.chunk=this.normalizeInput(e),this.file=n,this.type=this.constructor.type,this.globalOptions=this.options=t,this.localOptions=t[this.type],this.canTranslate=this.localOptions&&this.localOptions.translate}translate(){this.canTranslate&&(this.translated=this.translateBlock(this.raw,this.type))}get output(){return this.translated?this.translated:this.raw?Object.fromEntries(this.raw):void 0}translateBlock(e,t){let n=ie.get(t),r=j.get(t),s=D.get(t),a=this.options[t],o=a.reviveValues&&!!n,d=a.translateValues&&!!r,u=a.translateKeys&&!!s,c={};for(let[l,p]of e)o&&n.has(l)?p=n.get(l)(p):d&&r.has(l)&&(p=this.translateValue(p,r.get(l))),u&&s.has(l)&&(l=s.get(l)||l),c[l]=p;return c}translateValue(e,t){return t[e]||t.DEFAULT||e}assignToOutput(e,t){this.assignObjectToOutput(e,this.constructor.type,t)}assignObjectToOutput(e,t,n){if(this.globalOptions.mergeOutput)return Object.assign(e,n);e[t]?Object.assign(e[t],n):e[t]=n}}g(W,"headerLength",4),g(W,"type",void 0),g(W,"multiSegment",!1),g(W,"canHandle",()=>!1);function Zn(i){return i===192||i===194||i===196||i===219||i===221||i===218||i===254}function ei(i){return i>=224&&i<=239}function ti(i,e,t){for(let[n,r]of R)if(r.canHandle(i,e,t))return n}class Ot extends Xe{constructor(...e){super(...e),g(this,"appSegments",[]),g(this,"jpegSegments",[]),g(this,"unknownSegments",[])}static canHandle(e,t){return t===65496}async parse(){await this.findAppSegments(),await this.readSegments(this.appSegments),this.mergeMultiSegments(),this.createParsers(this.mergedAppSegments||this.appSegments)}setupSegmentFinderArgs(e){e===!0?(this.findAll=!0,this.wanted=new Set(R.keyList())):(e=e===void 0?R.keyList().filter(t=>this.options[t].enabled):e.filter(t=>this.options[t].enabled&&R.has(t)),this.findAll=!1,this.remaining=new Set(e),this.wanted=new Set(e)),this.unfinishedMultiSegment=!1}async findAppSegments(e=0,t){this.setupSegmentFinderArgs(t);let{file:n,findAll:r,wanted:s,remaining:a}=this;if(!r&&this.file.chunked&&(r=Array.from(s).some(o=>{let d=R.get(o),u=this.options[o];return d.multiSegment&&u.multiSegment}),r&&await this.file.readWhole()),e=this.findAppSegmentsInRange(e,n.byteLength),!this.options.onlyTiff&&n.chunked){let o=!1;for(;a.size>0&&!o&&(n.canReadNextChunk||this.unfinishedMultiSegment);){let{nextChunkOffset:d}=n,u=this.appSegments.some(c=>!this.file.available(c.offset||c.start,c.length||c.size));if(o=e>d&&!u?!await n.readNextChunk(e):!await n.readNextChunk(d),(e=this.findAppSegmentsInRange(e,n.byteLength))===void 0)return}}}findAppSegmentsInRange(e,t){t-=2;let n,r,s,a,o,d,{file:u,findAll:c,wanted:l,remaining:p,options:h}=this;for(;e<t;e++)if(u.getUint8(e)===255){if(n=u.getUint8(e+1),ei(n)){if(r=u.getUint16(e+2),s=ti(u,e,r),s&&l.has(s)&&(a=R.get(s),o=a.findPosition(u,e),d=h[s],o.type=s,this.appSegments.push(o),!c&&(a.multiSegment&&d.multiSegment?(this.unfinishedMultiSegment=o.chunkNumber<o.chunkCount,this.unfinishedMultiSegment||p.delete(s)):p.delete(s),p.size===0)))break;h.recordUnknownSegments&&(o=W.findPosition(u,e),o.marker=n,this.unknownSegments.push(o)),e+=r+1}else if(Zn(n)){if(r=u.getUint16(e+2),n===218&&h.stopAfterSos!==!1)return;h.recordJpegSegments&&this.jpegSegments.push({offset:e,length:r,marker:n}),e+=r+1}}return e}mergeMultiSegments(){if(!this.appSegments.some(t=>t.multiSegment))return;let e=function(t,n){let r,s,a,o=new Map;for(let d=0;d<t.length;d++)r=t[d],s=r[n],o.has(s)?a=o.get(s):o.set(s,a=[]),a.push(r);return Array.from(o)}(this.appSegments,"type");this.mergedAppSegments=e.map(([t,n])=>{let r=R.get(t,this.options);return r.handleMultiSegments?{type:t,chunk:r.handleMultiSegments(n)}:n[0]})}getSegment(e){return this.appSegments.find(t=>t.type===e)}async getOrFindSegment(e){let t=this.getSegment(e);return t===void 0&&(await this.findAppSegments(0,[e]),t=this.getSegment(e)),t}}g(Ot,"type","jpeg"),Y.set("jpeg",Ot);const ni=[void 0,1,1,2,4,8,1,1,2,4,8,4,8,4];class ii extends W{parseHeader(){var e=this.chunk.getUint16();e===18761?this.le=!0:e===19789&&(this.le=!1),this.chunk.le=this.le,this.headerParsed=!0}parseTags(e,t,n=new Map){let{pick:r,skip:s}=this.options[t];r=new Set(r);let a=r.size>0,o=s.size===0,d=this.chunk.getUint16(e);e+=2;for(let u=0;u<d;u++){let c=this.chunk.getUint16(e);if(a){if(r.has(c)&&(n.set(c,this.parseTag(e,c,t)),r.delete(c),r.size===0))break}else!o&&s.has(c)||n.set(c,this.parseTag(e,c,t));e+=12}return n}parseTag(e,t,n){let{chunk:r}=this,s=r.getUint16(e+2),a=r.getUint32(e+4),o=ni[s];if(o*a<=4?e+=8:e=r.getUint32(e+8),(s<1||s>13)&&O(`Invalid TIFF value type. block: ${n.toUpperCase()}, tag: ${t.toString(16)}, type: ${s}, offset ${e}`),e>r.byteLength&&O(`Invalid TIFF value offset. block: ${n.toUpperCase()}, tag: ${t.toString(16)}, type: ${s}, offset ${e} is outside of chunk size ${r.byteLength}`),s===1)return r.getUint8Array(e,a);if(s===2)return ce(r.getString(e,a));if(s===7)return r.getUint8Array(e,a);if(a===1)return this.parseTagValue(s,e);{let d=new(function(c){switch(c){case 1:return Uint8Array;case 3:return Uint16Array;case 4:return Uint32Array;case 5:return Array;case 6:return Int8Array;case 8:return Int16Array;case 9:return Int32Array;case 10:return Array;case 11:return Float32Array;case 12:return Float64Array;default:return Array}}(s))(a),u=o;for(let c=0;c<a;c++)d[c]=this.parseTagValue(s,e),e+=u;return d}}parseTagValue(e,t){let{chunk:n}=this;switch(e){case 1:return n.getUint8(t);case 3:return n.getUint16(t);case 4:return n.getUint32(t);case 5:return n.getUint32(t)/n.getUint32(t+4);case 6:return n.getInt8(t);case 8:return n.getInt16(t);case 9:return n.getInt32(t);case 10:return n.getInt32(t)/n.getInt32(t+4);case 11:return n.getFloat(t);case 12:return n.getDouble(t);case 13:return n.getUint32(t);default:O(`Invalid tiff type ${e}`)}}}class nt extends ii{static canHandle(e,t){return e.getUint8(t+1)===225&&e.getUint32(t+4)===1165519206&&e.getUint16(t+8)===0}async parse(){this.parseHeader();let{options:e}=this;return e.ifd0.enabled&&await this.parseIfd0Block(),e.exif.enabled&&await this.safeParse("parseExifBlock"),e.gps.enabled&&await this.safeParse("parseGpsBlock"),e.interop.enabled&&await this.safeParse("parseInteropBlock"),e.ifd1.enabled&&await this.safeParse("parseThumbnailBlock"),this.createOutput()}safeParse(e){let t=this[e]();return t.catch!==void 0&&(t=t.catch(this.handleError)),t}findIfd0Offset(){this.ifd0Offset===void 0&&(this.ifd0Offset=this.chunk.getUint32(4))}findIfd1Offset(){if(this.ifd1Offset===void 0){this.findIfd0Offset();let e=this.chunk.getUint16(this.ifd0Offset),t=this.ifd0Offset+2+12*e;this.ifd1Offset=this.chunk.getUint32(t)}}parseBlock(e,t){let n=new Map;return this[t]=n,this.parseTags(e,t,n),n}async parseIfd0Block(){if(this.ifd0)return;let{file:e}=this;this.findIfd0Offset(),this.ifd0Offset<8&&O("Malformed EXIF data"),!e.chunked&&this.ifd0Offset>e.byteLength&&O(`IFD0 offset points to outside of file.
this.ifd0Offset: ${this.ifd0Offset}, file.byteLength: ${e.byteLength}`),e.tiff&&await e.ensureChunk(this.ifd0Offset,dt(this.options));let t=this.parseBlock(this.ifd0Offset,"ifd0");return t.size!==0?(this.exifOffset=t.get(34665),this.interopOffset=t.get(40965),this.gpsOffset=t.get(34853),this.xmp=t.get(700),this.iptc=t.get(33723),this.icc=t.get(34675),this.options.sanitize&&(t.delete(34665),t.delete(40965),t.delete(34853),t.delete(700),t.delete(33723),t.delete(34675)),t):void 0}async parseExifBlock(){if(this.exif||(this.ifd0||await this.parseIfd0Block(),this.exifOffset===void 0))return;this.file.tiff&&await this.file.ensureChunk(this.exifOffset,dt(this.options));let e=this.parseBlock(this.exifOffset,"exif");return this.interopOffset||(this.interopOffset=e.get(40965)),this.makerNote=e.get(37500),this.userComment=e.get(37510),this.options.sanitize&&(e.delete(40965),e.delete(37500),e.delete(37510)),this.unpack(e,41728),this.unpack(e,41729),e}unpack(e,t){let n=e.get(t);n&&n.length===1&&e.set(t,n[0])}async parseGpsBlock(){if(this.gps||(this.ifd0||await this.parseIfd0Block(),this.gpsOffset===void 0))return;let e=this.parseBlock(this.gpsOffset,"gps");return e&&e.has(2)&&e.has(4)&&(e.set("latitude",Mt(...e.get(2),e.get(1))),e.set("longitude",Mt(...e.get(4),e.get(3)))),e}async parseInteropBlock(){if(!this.interop&&(this.ifd0||await this.parseIfd0Block(),this.interopOffset!==void 0||this.exif||await this.parseExifBlock(),this.interopOffset!==void 0))return this.parseBlock(this.interopOffset,"interop")}async parseThumbnailBlock(e=!1){if(!this.ifd1&&!this.ifd1Parsed&&(!this.options.mergeOutput||e))return this.findIfd1Offset(),this.ifd1Offset>0&&(this.parseBlock(this.ifd1Offset,"ifd1"),this.ifd1Parsed=!0),this.ifd1}async extractThumbnail(){if(this.headerParsed||this.parseHeader(),this.ifd1Parsed||await this.parseThumbnailBlock(!0),this.ifd1===void 0)return;let e=this.ifd1.get(513),t=this.ifd1.get(514);return this.chunk.getUint8Array(e,t)}get image(){return this.ifd0}get thumbnail(){return this.ifd1}createOutput(){let e,t,n,r={};for(t of v)if(e=this[t],!Qt(e))if(n=this.canTranslate?this.translateBlock(e,t):Object.fromEntries(e),this.options.mergeOutput){if(t==="ifd1")continue;Object.assign(r,n)}else r[t]=n;return this.makerNote&&(r.makerNote=this.makerNote),this.userComment&&(r.userComment=this.userComment),r}assignToOutput(e,t){if(this.globalOptions.mergeOutput)Object.assign(e,t);else for(let[n,r]of Object.entries(t))this.assignObjectToOutput(e,n,r)}}function Mt(i,e,t,n){var r=i+e/60+t/3600;return n!=="S"&&n!=="W"||(r*=-1),r}g(nt,"type","tiff"),g(nt,"headerLength",10),R.set("tiff",nt);var ri=Object.freeze({__proto__:null,default:Qn,Exifr:re,fileParsers:Y,segmentParsers:R,fileReaders:$,tagKeys:D,tagValues:j,tagRevivers:ie,createDictionary:I,extendDictionary:ge,fetchUrlAsArrayBuffer:pe,readBlobAsArrayBuffer:fe,chunkedProps:Z,otherSegments:be,segments:me,tiffBlocks:v,segmentsAndBlocks:ee,tiffExtractables:te,inheritables:xe,allFormatters:ne,Options:_e,parse:je});const mt={ifd0:!1,ifd1:!1,exif:!1,gps:!1,interop:!1,sanitize:!1,reviveValues:!0,translateKeys:!1,translateValues:!1,mergeOutput:!1},_t=Object.assign({},mt,{firstChunkSize:4e4,gps:[1,2,3,4]});async function tn(i){let e=new re(_t);await e.read(i);let t=await e.parse();if(t&&t.gps){let{latitude:n,longitude:r}=t.gps;return{latitude:n,longitude:r}}}const yt=Object.assign({},mt,{tiff:!1,ifd1:!0,mergeOutput:!1});async function nn(i){let e=new re(yt);await e.read(i);let t=await e.extractThumbnail();return t&&We?Ge.from(t):t}async function rn(i){let e=await this.thumbnail(i);if(e!==void 0){let t=new Blob([e]);return URL.createObjectURL(t)}}const wt=Object.assign({},mt,{firstChunkSize:4e4,ifd0:[274]});async function bt(i){let e=new re(wt);await e.read(i);let t=await e.parse();if(t&&t.ifd0)return t.ifd0[274]}const xt=Object.freeze({1:{dimensionSwapped:!1,scaleX:1,scaleY:1,deg:0,rad:0},2:{dimensionSwapped:!1,scaleX:-1,scaleY:1,deg:0,rad:0},3:{dimensionSwapped:!1,scaleX:1,scaleY:1,deg:180,rad:180*Math.PI/180},4:{dimensionSwapped:!1,scaleX:-1,scaleY:1,deg:180,rad:180*Math.PI/180},5:{dimensionSwapped:!0,scaleX:1,scaleY:-1,deg:90,rad:90*Math.PI/180},6:{dimensionSwapped:!0,scaleX:1,scaleY:1,deg:90,rad:90*Math.PI/180},7:{dimensionSwapped:!0,scaleX:1,scaleY:-1,deg:270,rad:270*Math.PI/180},8:{dimensionSwapped:!0,scaleX:1,scaleY:1,deg:270,rad:270*Math.PI/180}});let he=!0,ue=!0;if(typeof navigator=="object"){let i=navigator.userAgent;if(i.includes("iPad")||i.includes("iPhone")){let e=i.match(/OS (\d+)_(\d+)/);if(e){let[,t,n]=e;he=Number(t)+.1*Number(n)<13.4,ue=!1}}else if(i.includes("OS X 10")){let[,e]=i.match(/OS X 10[_.](\d+)/);he=ue=Number(e)<15}if(i.includes("Chrome/")){let[,e]=i.match(/Chrome\/(\d+)/);he=ue=Number(e)<81}else if(i.includes("Firefox/")){let[,e]=i.match(/Firefox\/(\d+)/);he=ue=Number(e)<77}}async function sn(i){let e=await bt(i);return Object.assign({canvas:he,css:ue},xt[e])}class si extends B{constructor(...e){super(...e),g(this,"ranges",new ai),this.byteLength!==0&&this.ranges.add(0,this.byteLength)}_tryExtend(e,t,n){if(e===0&&this.byteLength===0&&n){let r=new DataView(n.buffer||n,n.byteOffset,n.byteLength);this._swapDataView(r)}else{let r=e+t;if(r>this.byteLength){let{dataView:s}=this._extend(r);this._swapDataView(s)}}}_extend(e){let t;t=We?Ge.allocUnsafe(e):new Uint8Array(e);let n=new DataView(t.buffer,t.byteOffset,t.byteLength);return t.set(new Uint8Array(this.buffer,this.byteOffset,this.byteLength),0),{uintView:t,dataView:n}}subarray(e,t,n=!1){return t=t||this._lengthToEnd(e),n&&this._tryExtend(e,t),this.ranges.add(e,t),super.subarray(e,t)}set(e,t,n=!1){n&&this._tryExtend(t,e.byteLength,e);let r=super.set(e,t);return this.ranges.add(t,r.byteLength),r}async ensureChunk(e,t){this.chunked&&(this.ranges.available(e,t)||await this.readChunk(e,t))}available(e,t){return this.ranges.available(e,t)}}class ai{constructor(){g(this,"list",[])}get length(){return this.list.length}add(e,t,n=0){let r=e+t,s=this.list.filter(a=>Dt(e,a.offset,r)||Dt(e,a.end,r));if(s.length>0){e=Math.min(e,...s.map(o=>o.offset)),r=Math.max(r,...s.map(o=>o.end)),t=r-e;let a=s.shift();a.offset=e,a.length=t,a.end=r,this.list=this.list.filter(o=>!s.includes(o))}else this.list.push({offset:e,length:t,end:r})}available(e,t){let n=e+t;return this.list.some(r=>r.offset<=e&&n<=r.end)}}function Dt(i,e,t){return i<=e&&e<=t}class Ke extends si{constructor(e,t){super(0),g(this,"chunksRead",0),this.input=e,this.options=t}async readWhole(){this.chunked=!1,await this.readChunk(this.nextChunkOffset)}async readChunked(){this.chunked=!0,await this.readChunk(0,this.options.firstChunkSize)}async readNextChunk(e=this.nextChunkOffset){if(this.fullyRead)return this.chunksRead++,!1;let t=this.options.chunkSize,n=await this.readChunk(e,t);return!!n&&n.byteLength===t}async readChunk(e,t){if(this.chunksRead++,(t=this.safeWrapAddress(e,t))!==0)return this._readChunk(e,t)}safeWrapAddress(e,t){return this.size!==void 0&&e+t>this.size?Math.max(0,this.size-e):t}get nextChunkOffset(){if(this.ranges.list.length!==0)return this.ranges.list[0].length}get canReadNextChunk(){return this.chunksRead<this.options.chunkLimit}get fullyRead(){return this.size!==void 0&&this.nextChunkOffset===this.size}read(){return this.options.chunked?this.readChunked():this.readWhole()}close(){}}$.set("blob",class extends Ke{async readWhole(){this.chunked=!1;let i=await fe(this.input);this._swapArrayBuffer(i)}readChunked(){return this.chunked=!0,this.size=this.input.size,super.readChunked()}async _readChunk(i,e){let t=e?i+e:void 0,n=this.input.slice(i,t),r=await fe(n);return this.set(r,i,!0)}});var oi=Object.freeze({__proto__:null,default:ri,Exifr:re,fileParsers:Y,segmentParsers:R,fileReaders:$,tagKeys:D,tagValues:j,tagRevivers:ie,createDictionary:I,extendDictionary:ge,fetchUrlAsArrayBuffer:pe,readBlobAsArrayBuffer:fe,chunkedProps:Z,otherSegments:be,segments:me,tiffBlocks:v,segmentsAndBlocks:ee,tiffExtractables:te,inheritables:xe,allFormatters:ne,Options:_e,parse:je,gpsOnlyOptions:_t,gps:tn,thumbnailOnlyOptions:yt,thumbnail:nn,thumbnailUrl:rn,orientationOnlyOptions:wt,orientation:bt,rotations:xt,get rotateCanvas(){return he},get rotateCss(){return ue},rotation:sn});$.set("url",class extends Ke{async readWhole(){this.chunked=!1;let i=await pe(this.input);i instanceof ArrayBuffer?this._swapArrayBuffer(i):i instanceof Uint8Array&&this._swapBuffer(i)}async _readChunk(i,e){let t=e?i+e-1:void 0,n=this.options.httpHeaders||{};(i||t)&&(n.range=`bytes=${[i,t].join("-")}`);let r=await ft(this.input,{headers:n}),s=await r.arrayBuffer(),a=s.byteLength;if(r.status!==416)return a!==e&&(this.size=i+a),this.set(s,i,!0)}});B.prototype.getUint64=function(i){let e=this.getUint32(i),t=this.getUint32(i+4);return e<1048575?e<<32|t:typeof Le!==void 0?(console.warn("Using BigInt because of type 64uint but JS can only handle 53b numbers."),Le(e)<<Le(32)|Le(t)):void O("Trying to read 64b value but JS can only handle 53b numbers.")};class li extends Xe{parseBoxes(e=0){let t=[];for(;e<this.file.byteLength-4;){let n=this.parseBoxHead(e);if(t.push(n),n.length===0)break;e+=n.length}return t}parseSubBoxes(e){e.boxes=this.parseBoxes(e.start)}findBox(e,t){return e.boxes===void 0&&this.parseSubBoxes(e),e.boxes.find(n=>n.kind===t)}parseBoxHead(e){let t=this.file.getUint32(e),n=this.file.getString(e+4,4),r=e+8;return t===1&&(t=this.file.getUint64(e+8),r+=8),{offset:e,length:t,kind:n,start:r}}parseBoxFullHead(e){if(e.version!==void 0)return;let t=this.file.getUint32(e.start);e.version=t>>24,e.start+=4}}class an extends li{static canHandle(e,t){if(t!==0)return!1;let n=e.getUint16(2);if(n>50)return!1;let r=16,s=[];for(;r<n;)s.push(e.getString(r,4)),r+=4;return s.includes(this.type)}async parse(){let e=this.file.getUint32(0),t=this.parseBoxHead(e);for(;t.kind!=="meta";)e+=t.length,await this.file.ensureChunk(e,16),t=this.parseBoxHead(e);await this.file.ensureChunk(t.offset,t.length),this.parseBoxFullHead(t),this.parseSubBoxes(t),this.options.icc.enabled&&await this.findIcc(t),this.options.tiff.enabled&&await this.findExif(t)}async registerSegment(e,t,n){await this.file.ensureChunk(t,n);let r=this.file.subarray(t,n);this.createParser(e,r)}async findIcc(e){let t=this.findBox(e,"iprp");if(t===void 0)return;let n=this.findBox(t,"ipco");if(n===void 0)return;let r=this.findBox(n,"colr");r!==void 0&&await this.registerSegment("icc",r.offset+12,r.length)}async findExif(e){let t=this.findBox(e,"iinf");if(t===void 0)return;let n=this.findBox(e,"iloc");if(n===void 0)return;let r=this.findExifLocIdInIinf(t),s=this.findExtentInIloc(n,r);if(s===void 0)return;let[a,o]=s;await this.file.ensureChunk(a,o);let d=4+this.file.getUint32(a);a+=d,o-=d,await this.registerSegment("tiff",a,o)}findExifLocIdInIinf(e){this.parseBoxFullHead(e);let t,n,r,s,a=e.start,o=this.file.getUint16(a);for(a+=2;o--;){if(t=this.parseBoxHead(a),this.parseBoxFullHead(t),n=t.start,t.version>=2&&(r=t.version===3?4:2,s=this.file.getString(n+r+2,4),s==="Exif"))return this.file.getUintBytes(n,r);a+=t.length}}get8bits(e){let t=this.file.getUint8(e);return[t>>4,15&t]}findExtentInIloc(e,t){this.parseBoxFullHead(e);let n=e.start,[r,s]=this.get8bits(n++),[a,o]=this.get8bits(n++),d=e.version===2?4:2,u=e.version===1||e.version===2?2:0,c=o+r+s,l=e.version===2?4:2,p=this.file.getUintBytes(n,l);for(n+=l;p--;){let h=this.file.getUintBytes(n,d);n+=d+u+2+a;let m=this.file.getUint16(n);if(n+=2,h===t)return m>1&&console.warn(`ILOC box has more than one extent but we're only processing one
Please create an issue at https://github.com/MikeKovarik/exifr with this file`),[this.file.getUintBytes(n+o,r),this.file.getUintBytes(n+o+r,s)];n+=m*c}}}class on extends an{}g(on,"type","heic");class Ft extends an{}g(Ft,"type","avif"),Y.set("heic",on),Y.set("avif",Ft),I(D,["ifd0","ifd1"],[[256,"ImageWidth"],[257,"ImageHeight"],[258,"BitsPerSample"],[259,"Compression"],[262,"PhotometricInterpretation"],[270,"ImageDescription"],[271,"Make"],[272,"Model"],[273,"StripOffsets"],[274,"Orientation"],[277,"SamplesPerPixel"],[278,"RowsPerStrip"],[279,"StripByteCounts"],[282,"XResolution"],[283,"YResolution"],[284,"PlanarConfiguration"],[296,"ResolutionUnit"],[301,"TransferFunction"],[305,"Software"],[306,"ModifyDate"],[315,"Artist"],[316,"HostComputer"],[317,"Predictor"],[318,"WhitePoint"],[319,"PrimaryChromaticities"],[513,"ThumbnailOffset"],[514,"ThumbnailLength"],[529,"YCbCrCoefficients"],[530,"YCbCrSubSampling"],[531,"YCbCrPositioning"],[532,"ReferenceBlackWhite"],[700,"ApplicationNotes"],[33432,"Copyright"],[33723,"IPTC"],[34665,"ExifIFD"],[34675,"ICC"],[34853,"GpsIFD"],[330,"SubIFD"],[40965,"InteropIFD"],[40091,"XPTitle"],[40092,"XPComment"],[40093,"XPAuthor"],[40094,"XPKeywords"],[40095,"XPSubject"]]),I(D,"exif",[[33434,"ExposureTime"],[33437,"FNumber"],[34850,"ExposureProgram"],[34852,"SpectralSensitivity"],[34855,"ISO"],[34858,"TimeZoneOffset"],[34859,"SelfTimerMode"],[34864,"SensitivityType"],[34865,"StandardOutputSensitivity"],[34866,"RecommendedExposureIndex"],[34867,"ISOSpeed"],[34868,"ISOSpeedLatitudeyyy"],[34869,"ISOSpeedLatitudezzz"],[36864,"ExifVersion"],[36867,"DateTimeOriginal"],[36868,"CreateDate"],[36873,"GooglePlusUploadCode"],[36880,"OffsetTime"],[36881,"OffsetTimeOriginal"],[36882,"OffsetTimeDigitized"],[37121,"ComponentsConfiguration"],[37122,"CompressedBitsPerPixel"],[37377,"ShutterSpeedValue"],[37378,"ApertureValue"],[37379,"BrightnessValue"],[37380,"ExposureCompensation"],[37381,"MaxApertureValue"],[37382,"SubjectDistance"],[37383,"MeteringMode"],[37384,"LightSource"],[37385,"Flash"],[37386,"FocalLength"],[37393,"ImageNumber"],[37394,"SecurityClassification"],[37395,"ImageHistory"],[37396,"SubjectArea"],[37500,"MakerNote"],[37510,"UserComment"],[37520,"SubSecTime"],[37521,"SubSecTimeOriginal"],[37522,"SubSecTimeDigitized"],[37888,"AmbientTemperature"],[37889,"Humidity"],[37890,"Pressure"],[37891,"WaterDepth"],[37892,"Acceleration"],[37893,"CameraElevationAngle"],[40960,"FlashpixVersion"],[40961,"ColorSpace"],[40962,"ExifImageWidth"],[40963,"ExifImageHeight"],[40964,"RelatedSoundFile"],[41483,"FlashEnergy"],[41486,"FocalPlaneXResolution"],[41487,"FocalPlaneYResolution"],[41488,"FocalPlaneResolutionUnit"],[41492,"SubjectLocation"],[41493,"ExposureIndex"],[41495,"SensingMethod"],[41728,"FileSource"],[41729,"SceneType"],[41730,"CFAPattern"],[41985,"CustomRendered"],[41986,"ExposureMode"],[41987,"WhiteBalance"],[41988,"DigitalZoomRatio"],[41989,"FocalLengthIn35mmFormat"],[41990,"SceneCaptureType"],[41991,"GainControl"],[41992,"Contrast"],[41993,"Saturation"],[41994,"Sharpness"],[41996,"SubjectDistanceRange"],[42016,"ImageUniqueID"],[42032,"OwnerName"],[42033,"SerialNumber"],[42034,"LensInfo"],[42035,"LensMake"],[42036,"LensModel"],[42037,"LensSerialNumber"],[42080,"CompositeImage"],[42081,"CompositeImageCount"],[42082,"CompositeImageExposureTimes"],[42240,"Gamma"],[59932,"Padding"],[59933,"OffsetSchema"],[65e3,"OwnerName"],[65001,"SerialNumber"],[65002,"Lens"],[65100,"RawFile"],[65101,"Converter"],[65102,"WhiteBalance"],[65105,"Exposure"],[65106,"Shadows"],[65107,"Brightness"],[65108,"Contrast"],[65109,"Saturation"],[65110,"Sharpness"],[65111,"Smoothness"],[65112,"MoireFilter"],[40965,"InteropIFD"]]),I(D,"gps",[[0,"GPSVersionID"],[1,"GPSLatitudeRef"],[2,"GPSLatitude"],[3,"GPSLongitudeRef"],[4,"GPSLongitude"],[5,"GPSAltitudeRef"],[6,"GPSAltitude"],[7,"GPSTimeStamp"],[8,"GPSSatellites"],[9,"GPSStatus"],[10,"GPSMeasureMode"],[11,"GPSDOP"],[12,"GPSSpeedRef"],[13,"GPSSpeed"],[14,"GPSTrackRef"],[15,"GPSTrack"],[16,"GPSImgDirectionRef"],[17,"GPSImgDirection"],[18,"GPSMapDatum"],[19,"GPSDestLatitudeRef"],[20,"GPSDestLatitude"],[21,"GPSDestLongitudeRef"],[22,"GPSDestLongitude"],[23,"GPSDestBearingRef"],[24,"GPSDestBearing"],[25,"GPSDestDistanceRef"],[26,"GPSDestDistance"],[27,"GPSProcessingMethod"],[28,"GPSAreaInformation"],[29,"GPSDateStamp"],[30,"GPSDifferential"],[31,"GPSHPositioningError"]]),I(j,["ifd0","ifd1"],[[274,{1:"Horizontal (normal)",2:"Mirror horizontal",3:"Rotate 180",4:"Mirror vertical",5:"Mirror horizontal and rotate 270 CW",6:"Rotate 90 CW",7:"Mirror horizontal and rotate 90 CW",8:"Rotate 270 CW"}],[296,{1:"None",2:"inches",3:"cm"}]]);let Ae=I(j,"exif",[[34850,{0:"Not defined",1:"Manual",2:"Normal program",3:"Aperture priority",4:"Shutter priority",5:"Creative program",6:"Action program",7:"Portrait mode",8:"Landscape mode"}],[37121,{0:"-",1:"Y",2:"Cb",3:"Cr",4:"R",5:"G",6:"B"}],[37383,{0:"Unknown",1:"Average",2:"CenterWeightedAverage",3:"Spot",4:"MultiSpot",5:"Pattern",6:"Partial",255:"Other"}],[37384,{0:"Unknown",1:"Daylight",2:"Fluorescent",3:"Tungsten (incandescent light)",4:"Flash",9:"Fine weather",10:"Cloudy weather",11:"Shade",12:"Daylight fluorescent (D 5700 - 7100K)",13:"Day white fluorescent (N 4600 - 5400K)",14:"Cool white fluorescent (W 3900 - 4500K)",15:"White fluorescent (WW 3200 - 3700K)",17:"Standard light A",18:"Standard light B",19:"Standard light C",20:"D55",21:"D65",22:"D75",23:"D50",24:"ISO studio tungsten",255:"Other"}],[37385,{0:"Flash did not fire",1:"Flash fired",5:"Strobe return light not detected",7:"Strobe return light detected",9:"Flash fired, compulsory flash mode",13:"Flash fired, compulsory flash mode, return light not detected",15:"Flash fired, compulsory flash mode, return light detected",16:"Flash did not fire, compulsory flash mode",24:"Flash did not fire, auto mode",25:"Flash fired, auto mode",29:"Flash fired, auto mode, return light not detected",31:"Flash fired, auto mode, return light detected",32:"No flash function",65:"Flash fired, red-eye reduction mode",69:"Flash fired, red-eye reduction mode, return light not detected",71:"Flash fired, red-eye reduction mode, return light detected",73:"Flash fired, compulsory flash mode, red-eye reduction mode",77:"Flash fired, compulsory flash mode, red-eye reduction mode, return light not detected",79:"Flash fired, compulsory flash mode, red-eye reduction mode, return light detected",89:"Flash fired, auto mode, red-eye reduction mode",93:"Flash fired, auto mode, return light not detected, red-eye reduction mode",95:"Flash fired, auto mode, return light detected, red-eye reduction mode"}],[41495,{1:"Not defined",2:"One-chip color area sensor",3:"Two-chip color area sensor",4:"Three-chip color area sensor",5:"Color sequential area sensor",7:"Trilinear sensor",8:"Color sequential linear sensor"}],[41728,{1:"Film Scanner",2:"Reflection Print Scanner",3:"Digital Camera"}],[41729,{1:"Directly photographed"}],[41985,{0:"Normal",1:"Custom",2:"HDR (no original saved)",3:"HDR (original saved)",4:"Original (for HDR)",6:"Panorama",7:"Portrait HDR",8:"Portrait"}],[41986,{0:"Auto",1:"Manual",2:"Auto bracket"}],[41987,{0:"Auto",1:"Manual"}],[41990,{0:"Standard",1:"Landscape",2:"Portrait",3:"Night",4:"Other"}],[41991,{0:"None",1:"Low gain up",2:"High gain up",3:"Low gain down",4:"High gain down"}],[41996,{0:"Unknown",1:"Macro",2:"Close",3:"Distant"}],[42080,{0:"Unknown",1:"Not a Composite Image",2:"General Composite Image",3:"Composite Image Captured While Shooting"}]]);const Et={1:"No absolute unit of measurement",2:"Inch",3:"Centimeter"};Ae.set(37392,Et),Ae.set(41488,Et);const it={0:"Normal",1:"Low",2:"High"};function Nt(i){return typeof i=="object"&&i.length!==void 0?i[0]:i}function Ut(i){let e=Array.from(i).slice(1);return e[1]>15&&(e=e.map(t=>String.fromCharCode(t))),e[2]!=="0"&&e[2]!==0||e.pop(),e.join(".")}function rt(i){if(typeof i=="string"){var[e,t,n,r,s,a]=i.trim().split(/[-: ]/g).map(Number),o=new Date(e,t-1,n);return Number.isNaN(r)||Number.isNaN(s)||Number.isNaN(a)||(o.setHours(r),o.setMinutes(s),o.setSeconds(a)),Number.isNaN(+o)?i:o}}function Ce(i){if(typeof i=="string")return i;let e=[];if(i[1]===0&&i[i.length-1]===0)for(let t=0;t<i.length;t+=2)e.push(Bt(i[t+1],i[t]));else for(let t=0;t<i.length;t+=2)e.push(Bt(i[t],i[t+1]));return ce(String.fromCodePoint(...e))}function Bt(i,e){return i<<8|e}Ae.set(41992,it),Ae.set(41993,it),Ae.set(41994,it),I(ie,["ifd0","ifd1"],[[50827,function(i){return typeof i!="string"?Zt(i):i}],[306,rt],[40091,Ce],[40092,Ce],[40093,Ce],[40094,Ce],[40095,Ce]]),I(ie,"exif",[[40960,Ut],[36864,Ut],[36867,rt],[36868,rt],[40962,Nt],[40963,Nt]]),I(ie,"gps",[[0,i=>Array.from(i).join(".")],[7,i=>Array.from(i).join(":")]]);class st extends W{static canHandle(e,t){return e.getUint8(t+1)===225&&e.getUint32(t+4)===1752462448&&e.getString(t+4,20)==="http://ns.adobe.com/"}static headerLength(e,t){return e.getString(t+4,34)==="http://ns.adobe.com/xmp/extension/"?79:33}static findPosition(e,t){let n=super.findPosition(e,t);return n.multiSegment=n.extended=n.headerLength===79,n.multiSegment?(n.chunkCount=e.getUint8(t+72),n.chunkNumber=e.getUint8(t+76),e.getUint8(t+77)!==0&&n.chunkNumber++):(n.chunkCount=1/0,n.chunkNumber=-1),n}static handleMultiSegments(e){return e.map(t=>t.chunk.getString()).join("")}normalizeInput(e){return typeof e=="string"?e:B.from(e).getString()}parse(e=this.chunk){if(!this.localOptions.parse)return e;e=function(s){let a={},o={};for(let d of hn)a[d]=[],o[d]=0;return s.replace(ui,(d,u,c)=>{if(u==="<"){let l=++o[c];return a[c].push(l),`${d}#${l}`}return`${d}#${a[c].pop()}`})}(e);let t=we.findAll(e,"rdf","Description");t.length===0&&t.push(new we("rdf","Description",void 0,e));let n,r={};for(let s of t)for(let a of s.properties)n=hi(a.ns,r),ln(a,n);return function(s){let a;for(let o in s)a=s[o]=Ue(s[o]),a===void 0&&delete s[o];return Ue(s)}(r)}assignToOutput(e,t){if(this.localOptions.parse)for(let[n,r]of Object.entries(t))switch(n){case"tiff":this.assignObjectToOutput(e,"ifd0",r);break;case"exif":this.assignObjectToOutput(e,"exif",r);break;case"xmlns":break;default:this.assignObjectToOutput(e,n,r)}else e.xmp=t}}g(st,"type","xmp"),g(st,"multiSegment",!0),R.set("xmp",st);class Ve{static findAll(e){return dn(e,/([a-zA-Z0-9-]+):([a-zA-Z0-9-]+)=("[^"]*"|'[^']*')/gm).map(Ve.unpackMatch)}static unpackMatch(e){let t=e[1],n=e[2],r=e[3].slice(1,-1);return r=cn(r),new Ve(t,n,r)}constructor(e,t,n){this.ns=e,this.name=t,this.value=n}serialize(){return this.value}}class we{static findAll(e,t,n){if(t!==void 0||n!==void 0){t=t||"[\\w\\d-]+",n=n||"[\\w\\d-]+";var r=new RegExp(`<(${t}):(${n})(#\\d+)?((\\s+?[\\w\\d-:]+=("[^"]*"|'[^']*'))*\\s*)(\\/>|>([\\s\\S]*?)<\\/\\1:\\2\\3>)`,"gm")}else r=/<([\w\d-]+):([\w\d-]+)(#\d+)?((\s+?[\w\d-:]+=("[^"]*"|'[^']*'))*\s*)(\/>|>([\s\S]*?)<\/\1:\2\3>)/gm;return dn(e,r).map(we.unpackMatch)}static unpackMatch(e){let t=e[1],n=e[2],r=e[4],s=e[8];return new we(t,n,r,s)}constructor(e,t,n,r){this.ns=e,this.name=t,this.attrString=n,this.innerXml=r,this.attrs=Ve.findAll(n),this.children=we.findAll(r),this.value=this.children.length===0?cn(r):void 0,this.properties=[...this.attrs,...this.children]}get isPrimitive(){return this.value!==void 0&&this.attrs.length===0&&this.children.length===0}get isListContainer(){return this.children.length===1&&this.children[0].isList}get isList(){let{ns:e,name:t}=this;return e==="rdf"&&(t==="Seq"||t==="Bag"||t==="Alt")}get isListItem(){return this.ns==="rdf"&&this.name==="li"}serialize(){if(this.properties.length===0&&this.value===void 0)return;if(this.isPrimitive)return this.value;if(this.isListContainer)return this.children[0].serialize();if(this.isList)return ci(this.children.map(di));if(this.isListItem&&this.children.length===1&&this.attrs.length===0)return this.children[0].serialize();let e={};for(let t of this.properties)ln(t,e);return this.value!==void 0&&(e.value=this.value),Ue(e)}}function ln(i,e){let t=i.serialize();t!==void 0&&(e[i.name]=t)}var di=i=>i.serialize(),ci=i=>i.length===1?i[0]:i,hi=(i,e)=>e[i]?e[i]:e[i]={};function dn(i,e){let t,n=[];if(!i)return n;for(;(t=e.exec(i))!==null;)n.push(t);return n}function cn(i){if(function(n){return n==null||n==="null"||n==="undefined"||n===""||n.trim()===""}(i))return;let e=Number(i);if(!Number.isNaN(e))return e;let t=i.toLowerCase();return t==="true"||t!=="false"&&i.trim()}const hn=["rdf:li","rdf:Seq","rdf:Bag","rdf:Alt","rdf:Description"],ui=new RegExp(`(<|\\/)(${hn.join("|")})`,"g");var pi=Object.freeze({__proto__:null,default:oi,Exifr:re,fileParsers:Y,segmentParsers:R,fileReaders:$,tagKeys:D,tagValues:j,tagRevivers:ie,createDictionary:I,extendDictionary:ge,fetchUrlAsArrayBuffer:pe,readBlobAsArrayBuffer:fe,chunkedProps:Z,otherSegments:be,segments:me,tiffBlocks:v,segmentsAndBlocks:ee,tiffExtractables:te,inheritables:xe,allFormatters:ne,Options:_e,parse:je,gpsOnlyOptions:_t,gps:tn,thumbnailOnlyOptions:yt,thumbnail:nn,thumbnailUrl:rn,orientationOnlyOptions:wt,orientation:bt,rotations:xt,get rotateCanvas(){return he},get rotateCss(){return ue},rotation:sn});let Vt=Ne("fs",i=>i.promises);$.set("fs",class extends Ke{async readWhole(){this.chunked=!1,this.fs=await Vt;let i=await this.fs.readFile(this.input);this._swapBuffer(i)}async readChunked(){this.chunked=!0,this.fs=await Vt,await this.open(),await this.readChunk(0,this.options.firstChunkSize)}async open(){this.fh===void 0&&(this.fh=await this.fs.open(this.input,"r"),this.size=(await this.fh.stat(this.input)).size)}async _readChunk(i,e){this.fh===void 0&&await this.open(),i+e>this.size&&(e=this.size-i);var t=this.subarray(i,e,!0);return await this.fh.read(t.dataView,0,e,i),t}async close(){if(this.fh){let i=this.fh;this.fh=void 0,await i.close()}}});$.set("base64",class extends Ke{constructor(...i){super(...i),this.input=this.input.replace(/^data:([^;]+);base64,/gim,""),this.size=this.input.length/4*3,this.input.endsWith("==")?this.size-=2:this.input.endsWith("=")&&(this.size-=1)}async _readChunk(i,e){let t,n,r=this.input;i===void 0?(i=0,t=0,n=0):(t=4*Math.floor(i/3),n=i-t/4*3),e===void 0&&(e=this.size);let s=i+e,a=t+4*Math.ceil(s/3);r=r.slice(t,a);let o=Math.min(e,this.size-i);if(We){let d=Ge.from(r,"base64").slice(n,n+o);return this.set(d,i,!0)}{let d=this.subarray(i,o,!0),u=atob(r),c=d.toUint8();for(let l=0;l<o;l++)c[l]=u.charCodeAt(n+l);return d}}});class zt extends Xe{static canHandle(e,t){return t===18761||t===19789}extendOptions(e){let{ifd0:t,xmp:n,iptc:r,icc:s}=e;n.enabled&&t.deps.add(700),r.enabled&&t.deps.add(33723),s.enabled&&t.deps.add(34675),t.finalizeFilters()}async parse(){let{tiff:e,xmp:t,iptc:n,icc:r}=this.options;if(e.enabled||t.enabled||n.enabled||r.enabled){let s=Math.max(dt(this.options),this.options.chunkSize);await this.file.ensureChunk(0,s),this.createParser("tiff",this.file),this.parsers.tiff.parseHeader(),await this.parsers.tiff.parseIfd0Block(),this.adaptTiffPropAsSegment("xmp"),this.adaptTiffPropAsSegment("iptc"),this.adaptTiffPropAsSegment("icc")}}adaptTiffPropAsSegment(e){if(this.parsers.tiff[e]){let t=this.parsers.tiff[e];this.injectSegment(e,t)}}}g(zt,"type","tiff"),Y.set("tiff",zt);let fi=Ne("zlib");const gi=["ihdr","iccp","text","itxt","exif"];class Ht extends Xe{constructor(...e){super(...e),g(this,"catchError",t=>this.errors.push(t)),g(this,"metaChunks",[]),g(this,"unknownChunks",[])}static canHandle(e,t){return t===35152&&e.getUint32(0)===2303741511&&e.getUint32(4)===218765834}async parse(){let{file:e}=this;await this.findPngChunksInRange(8,e.byteLength),await this.readSegments(this.metaChunks),this.findIhdr(),this.parseTextChunks(),await this.findExif().catch(this.catchError),await this.findXmp().catch(this.catchError),await this.findIcc().catch(this.catchError)}async findPngChunksInRange(e,t){let{file:n}=this;for(;e<t;){let r=n.getUint32(e),s=n.getUint32(e+4),a=n.getString(e+4,4).toLowerCase(),o=r+4+4+4,d={type:a,offset:e,length:o,start:e+4+4,size:r,marker:s};gi.includes(a)?this.metaChunks.push(d):this.unknownChunks.push(d),e+=o}}parseTextChunks(){let e=this.metaChunks.filter(t=>t.type==="text");for(let t of e){let[n,r]=this.file.getString(t.start,t.size).split("\0");this.injectKeyValToIhdr(n,r)}}injectKeyValToIhdr(e,t){let n=this.parsers.ihdr;n&&n.raw.set(e,t)}findIhdr(){let e=this.metaChunks.find(t=>t.type==="ihdr");e&&this.options.ihdr.enabled!==!1&&this.createParser("ihdr",e.chunk)}async findExif(){let e=this.metaChunks.find(t=>t.type==="exif");e&&this.injectSegment("tiff",e.chunk)}async findXmp(){let e=this.metaChunks.filter(t=>t.type==="itxt");for(let t of e)t.chunk.getString(0,17)==="XML:com.adobe.xmp"&&this.injectSegment("xmp",t.chunk)}async findIcc(){let e=this.metaChunks.find(o=>o.type==="iccp");if(!e)return;let{chunk:t}=e,n=t.getUint8Array(0,81),r=0;for(;r<80&&n[r]!==0;)r++;let s=r+2,a=t.getString(0,r);if(this.injectKeyValToIhdr("ProfileName",a),Ee){let o=await fi,d=t.getUint8Array(s);d=o.inflateSync(d),this.injectSegment("icc",d)}}}g(Ht,"type","png"),Y.set("png",Ht),I(D,"interop",[[1,"InteropIndex"],[2,"InteropVersion"],[4096,"RelatedImageFileFormat"],[4097,"RelatedImageWidth"],[4098,"RelatedImageHeight"]]),ge(D,"ifd0",[[11,"ProcessingSoftware"],[254,"SubfileType"],[255,"OldSubfileType"],[263,"Thresholding"],[264,"CellWidth"],[265,"CellLength"],[266,"FillOrder"],[269,"DocumentName"],[280,"MinSampleValue"],[281,"MaxSampleValue"],[285,"PageName"],[286,"XPosition"],[287,"YPosition"],[290,"GrayResponseUnit"],[297,"PageNumber"],[321,"HalftoneHints"],[322,"TileWidth"],[323,"TileLength"],[332,"InkSet"],[337,"TargetPrinter"],[18246,"Rating"],[18249,"RatingPercent"],[33550,"PixelScale"],[34264,"ModelTransform"],[34377,"PhotoshopSettings"],[50706,"DNGVersion"],[50707,"DNGBackwardVersion"],[50708,"UniqueCameraModel"],[50709,"LocalizedCameraModel"],[50736,"DNGLensInfo"],[50739,"ShadowScale"],[50740,"DNGPrivateData"],[33920,"IntergraphMatrix"],[33922,"ModelTiePoint"],[34118,"SEMInfo"],[34735,"GeoTiffDirectory"],[34736,"GeoTiffDoubleParams"],[34737,"GeoTiffAsciiParams"],[50341,"PrintIM"],[50721,"ColorMatrix1"],[50722,"ColorMatrix2"],[50723,"CameraCalibration1"],[50724,"CameraCalibration2"],[50725,"ReductionMatrix1"],[50726,"ReductionMatrix2"],[50727,"AnalogBalance"],[50728,"AsShotNeutral"],[50729,"AsShotWhiteXY"],[50730,"BaselineExposure"],[50731,"BaselineNoise"],[50732,"BaselineSharpness"],[50734,"LinearResponseLimit"],[50735,"CameraSerialNumber"],[50741,"MakerNoteSafety"],[50778,"CalibrationIlluminant1"],[50779,"CalibrationIlluminant2"],[50781,"RawDataUniqueID"],[50827,"OriginalRawFileName"],[50828,"OriginalRawFileData"],[50831,"AsShotICCProfile"],[50832,"AsShotPreProfileMatrix"],[50833,"CurrentICCProfile"],[50834,"CurrentPreProfileMatrix"],[50879,"ColorimetricReference"],[50885,"SRawType"],[50898,"PanasonicTitle"],[50899,"PanasonicTitle2"],[50931,"CameraCalibrationSig"],[50932,"ProfileCalibrationSig"],[50933,"ProfileIFD"],[50934,"AsShotProfileName"],[50936,"ProfileName"],[50937,"ProfileHueSatMapDims"],[50938,"ProfileHueSatMapData1"],[50939,"ProfileHueSatMapData2"],[50940,"ProfileToneCurve"],[50941,"ProfileEmbedPolicy"],[50942,"ProfileCopyright"],[50964,"ForwardMatrix1"],[50965,"ForwardMatrix2"],[50966,"PreviewApplicationName"],[50967,"PreviewApplicationVersion"],[50968,"PreviewSettingsName"],[50969,"PreviewSettingsDigest"],[50970,"PreviewColorSpace"],[50971,"PreviewDateTime"],[50972,"RawImageDigest"],[50973,"OriginalRawFileDigest"],[50981,"ProfileLookTableDims"],[50982,"ProfileLookTableData"],[51043,"TimeCodes"],[51044,"FrameRate"],[51058,"TStop"],[51081,"ReelName"],[51089,"OriginalDefaultFinalSize"],[51090,"OriginalBestQualitySize"],[51091,"OriginalDefaultCropSize"],[51105,"CameraLabel"],[51107,"ProfileHueSatMapEncoding"],[51108,"ProfileLookTableEncoding"],[51109,"BaselineExposureOffset"],[51110,"DefaultBlackRender"],[51111,"NewRawImageDigest"],[51112,"RawToPreviewGain"]]);let Gt=[[273,"StripOffsets"],[279,"StripByteCounts"],[288,"FreeOffsets"],[289,"FreeByteCounts"],[291,"GrayResponseCurve"],[292,"T4Options"],[293,"T6Options"],[300,"ColorResponseUnit"],[320,"ColorMap"],[324,"TileOffsets"],[325,"TileByteCounts"],[326,"BadFaxLines"],[327,"CleanFaxData"],[328,"ConsecutiveBadFaxLines"],[330,"SubIFD"],[333,"InkNames"],[334,"NumberofInks"],[336,"DotRange"],[338,"ExtraSamples"],[339,"SampleFormat"],[340,"SMinSampleValue"],[341,"SMaxSampleValue"],[342,"TransferRange"],[343,"ClipPath"],[344,"XClipPathUnits"],[345,"YClipPathUnits"],[346,"Indexed"],[347,"JPEGTables"],[351,"OPIProxy"],[400,"GlobalParametersIFD"],[401,"ProfileType"],[402,"FaxProfile"],[403,"CodingMethods"],[404,"VersionYear"],[405,"ModeNumber"],[433,"Decode"],[434,"DefaultImageColor"],[435,"T82Options"],[437,"JPEGTables"],[512,"JPEGProc"],[515,"JPEGRestartInterval"],[517,"JPEGLosslessPredictors"],[518,"JPEGPointTransforms"],[519,"JPEGQTables"],[520,"JPEGDCTables"],[521,"JPEGACTables"],[559,"StripRowCounts"],[999,"USPTOMiscellaneous"],[18247,"XP_DIP_XML"],[18248,"StitchInfo"],[28672,"SonyRawFileType"],[28688,"SonyToneCurve"],[28721,"VignettingCorrection"],[28722,"VignettingCorrParams"],[28724,"ChromaticAberrationCorrection"],[28725,"ChromaticAberrationCorrParams"],[28726,"DistortionCorrection"],[28727,"DistortionCorrParams"],[29895,"SonyCropTopLeft"],[29896,"SonyCropSize"],[32781,"ImageID"],[32931,"WangTag1"],[32932,"WangAnnotation"],[32933,"WangTag3"],[32934,"WangTag4"],[32953,"ImageReferencePoints"],[32954,"RegionXformTackPoint"],[32955,"WarpQuadrilateral"],[32956,"AffineTransformMat"],[32995,"Matteing"],[32996,"DataType"],[32997,"ImageDepth"],[32998,"TileDepth"],[33300,"ImageFullWidth"],[33301,"ImageFullHeight"],[33302,"TextureFormat"],[33303,"WrapModes"],[33304,"FovCot"],[33305,"MatrixWorldToScreen"],[33306,"MatrixWorldToCamera"],[33405,"Model2"],[33421,"CFARepeatPatternDim"],[33422,"CFAPattern2"],[33423,"BatteryLevel"],[33424,"KodakIFD"],[33445,"MDFileTag"],[33446,"MDScalePixel"],[33447,"MDColorTable"],[33448,"MDLabName"],[33449,"MDSampleInfo"],[33450,"MDPrepDate"],[33451,"MDPrepTime"],[33452,"MDFileUnits"],[33589,"AdventScale"],[33590,"AdventRevision"],[33628,"UIC1Tag"],[33629,"UIC2Tag"],[33630,"UIC3Tag"],[33631,"UIC4Tag"],[33918,"IntergraphPacketData"],[33919,"IntergraphFlagRegisters"],[33921,"INGRReserved"],[34016,"Site"],[34017,"ColorSequence"],[34018,"IT8Header"],[34019,"RasterPadding"],[34020,"BitsPerRunLength"],[34021,"BitsPerExtendedRunLength"],[34022,"ColorTable"],[34023,"ImageColorIndicator"],[34024,"BackgroundColorIndicator"],[34025,"ImageColorValue"],[34026,"BackgroundColorValue"],[34027,"PixelIntensityRange"],[34028,"TransparencyIndicator"],[34029,"ColorCharacterization"],[34030,"HCUsage"],[34031,"TrapIndicator"],[34032,"CMYKEquivalent"],[34152,"AFCP_IPTC"],[34232,"PixelMagicJBIGOptions"],[34263,"JPLCartoIFD"],[34306,"WB_GRGBLevels"],[34310,"LeafData"],[34687,"TIFF_FXExtensions"],[34688,"MultiProfiles"],[34689,"SharedData"],[34690,"T88Options"],[34732,"ImageLayer"],[34750,"JBIGOptions"],[34856,"Opto-ElectricConvFactor"],[34857,"Interlace"],[34908,"FaxRecvParams"],[34909,"FaxSubAddress"],[34910,"FaxRecvTime"],[34929,"FedexEDR"],[34954,"LeafSubIFD"],[37387,"FlashEnergy"],[37388,"SpatialFrequencyResponse"],[37389,"Noise"],[37390,"FocalPlaneXResolution"],[37391,"FocalPlaneYResolution"],[37392,"FocalPlaneResolutionUnit"],[37397,"ExposureIndex"],[37398,"TIFF-EPStandardID"],[37399,"SensingMethod"],[37434,"CIP3DataFile"],[37435,"CIP3Sheet"],[37436,"CIP3Side"],[37439,"StoNits"],[37679,"MSDocumentText"],[37680,"MSPropertySetStorage"],[37681,"MSDocumentTextPosition"],[37724,"ImageSourceData"],[40965,"InteropIFD"],[40976,"SamsungRawPointersOffset"],[40977,"SamsungRawPointersLength"],[41217,"SamsungRawByteOrder"],[41218,"SamsungRawUnknown"],[41484,"SpatialFrequencyResponse"],[41485,"Noise"],[41489,"ImageNumber"],[41490,"SecurityClassification"],[41491,"ImageHistory"],[41494,"TIFF-EPStandardID"],[41995,"DeviceSettingDescription"],[42112,"GDALMetadata"],[42113,"GDALNoData"],[44992,"ExpandSoftware"],[44993,"ExpandLens"],[44994,"ExpandFilm"],[44995,"ExpandFilterLens"],[44996,"ExpandScanner"],[44997,"ExpandFlashLamp"],[46275,"HasselbladRawImage"],[48129,"PixelFormat"],[48130,"Transformation"],[48131,"Uncompressed"],[48132,"ImageType"],[48256,"ImageWidth"],[48257,"ImageHeight"],[48258,"WidthResolution"],[48259,"HeightResolution"],[48320,"ImageOffset"],[48321,"ImageByteCount"],[48322,"AlphaOffset"],[48323,"AlphaByteCount"],[48324,"ImageDataDiscard"],[48325,"AlphaDataDiscard"],[50215,"OceScanjobDesc"],[50216,"OceApplicationSelector"],[50217,"OceIDNumber"],[50218,"OceImageLogic"],[50255,"Annotations"],[50459,"HasselbladExif"],[50547,"OriginalFileName"],[50560,"USPTOOriginalContentType"],[50656,"CR2CFAPattern"],[50710,"CFAPlaneColor"],[50711,"CFALayout"],[50712,"LinearizationTable"],[50713,"BlackLevelRepeatDim"],[50714,"BlackLevel"],[50715,"BlackLevelDeltaH"],[50716,"BlackLevelDeltaV"],[50717,"WhiteLevel"],[50718,"DefaultScale"],[50719,"DefaultCropOrigin"],[50720,"DefaultCropSize"],[50733,"BayerGreenSplit"],[50737,"ChromaBlurRadius"],[50738,"AntiAliasStrength"],[50752,"RawImageSegmentation"],[50780,"BestQualityScale"],[50784,"AliasLayerMetadata"],[50829,"ActiveArea"],[50830,"MaskedAreas"],[50935,"NoiseReductionApplied"],[50974,"SubTileBlockSize"],[50975,"RowInterleaveFactor"],[51008,"OpcodeList1"],[51009,"OpcodeList2"],[51022,"OpcodeList3"],[51041,"NoiseProfile"],[51114,"CacheVersion"],[51125,"DefaultUserCrop"],[51157,"NikonNEFInfo"],[65024,"KdcIFD"]];ge(D,"ifd0",Gt),ge(D,"exif",Gt),I(j,"gps",[[23,{M:"Magnetic North",T:"True North"}],[25,{K:"Kilometers",M:"Miles",N:"Nautical Miles"}]]);class at extends W{static canHandle(e,t){return e.getUint8(t+1)===224&&e.getUint32(t+4)===1246120262&&e.getUint8(t+8)===0}parse(){return this.parseTags(),this.translate(),this.output}parseTags(){this.raw=new Map([[0,this.chunk.getUint16(0)],[2,this.chunk.getUint8(2)],[3,this.chunk.getUint16(3)],[5,this.chunk.getUint16(5)],[7,this.chunk.getUint8(7)],[8,this.chunk.getUint8(8)]])}}g(at,"type","jfif"),g(at,"headerLength",9),R.set("jfif",at),I(D,"jfif",[[0,"JFIFVersion"],[2,"ResolutionUnit"],[3,"XResolution"],[5,"YResolution"],[7,"ThumbnailWidth"],[8,"ThumbnailHeight"]]);class Wt extends W{parse(){return this.parseTags(),this.translate(),this.output}parseTags(){this.raw=new Map([[0,this.chunk.getUint32(0)],[4,this.chunk.getUint32(4)],[8,this.chunk.getUint8(8)],[9,this.chunk.getUint8(9)],[10,this.chunk.getUint8(10)],[11,this.chunk.getUint8(11)],[12,this.chunk.getUint8(12)],...Array.from(this.raw)])}}g(Wt,"type","ihdr"),R.set("ihdr",Wt),I(D,"ihdr",[[0,"ImageWidth"],[4,"ImageHeight"],[8,"BitDepth"],[9,"ColorType"],[10,"Compression"],[11,"Filter"],[12,"Interlace"]]),I(j,"ihdr",[[9,{0:"Grayscale",2:"RGB",3:"Palette",4:"Grayscale with Alpha",6:"RGB with Alpha",DEFAULT:"Unknown"}],[10,{0:"Deflate/Inflate",DEFAULT:"Unknown"}],[11,{0:"Adaptive",DEFAULT:"Unknown"}],[12,{0:"Noninterlaced",1:"Adam7 Interlace",DEFAULT:"Unknown"}]]);class Me extends W{static canHandle(e,t){return e.getUint8(t+1)===226&&e.getUint32(t+4)===1229144927}static findPosition(e,t){let n=super.findPosition(e,t);return n.chunkNumber=e.getUint8(t+16),n.chunkCount=e.getUint8(t+17),n.multiSegment=n.chunkCount>1,n}static handleMultiSegments(e){return function(t){let n=function(r){let s=r[0].constructor,a=0;for(let u of r)a+=u.length;let o=new s(a),d=0;for(let u of r)o.set(u,d),d+=u.length;return o}(t.map(r=>r.chunk.toUint8()));return new B(n)}(e)}parse(){return this.raw=new Map,this.parseHeader(),this.parseTags(),this.translate(),this.output}parseHeader(){let{raw:e}=this;this.chunk.byteLength<84&&O("ICC header is too short");for(let[t,n]of Object.entries(mi)){t=parseInt(t,10);let r=n(this.chunk,t);r!=="\0\0\0\0"&&e.set(t,r)}}parseTags(){let e,t,n,r,s,{raw:a}=this,o=this.chunk.getUint32(128),d=132,u=this.chunk.byteLength;for(;o--;){if(e=this.chunk.getString(d,4),t=this.chunk.getUint32(d+4),n=this.chunk.getUint32(d+8),r=this.chunk.getString(t,4),t+n>u)return void console.warn("reached the end of the first ICC chunk. Enable options.tiff.multiSegment to read all ICC segments.");s=this.parseTag(r,t,n),s!==void 0&&s!=="\0\0\0\0"&&a.set(e,s),d+=12}}parseTag(e,t,n){switch(e){case"desc":return this.parseDesc(t);case"mluc":return this.parseMluc(t);case"text":return this.parseText(t,n);case"sig ":return this.parseSig(t)}if(!(t+n>this.chunk.byteLength))return this.chunk.getUint8Array(t,n)}parseDesc(e){let t=this.chunk.getUint32(e+8)-1;return ce(this.chunk.getString(e+12,t))}parseText(e,t){return ce(this.chunk.getString(e+8,t-8))}parseSig(e){return ce(this.chunk.getString(e+8,4))}parseMluc(e){let{chunk:t}=this,n=t.getUint32(e+8),r=t.getUint32(e+12),s=e+16,a=[];for(let o=0;o<n;o++){let d=t.getString(s+0,2),u=t.getString(s+2,2),c=t.getUint32(s+4),l=t.getUint32(s+8)+e,p=ce(t.getUnicodeString(l,c));a.push({lang:d,country:u,text:p}),s+=r}return n===1?a[0].text:a}translateValue(e,t){return typeof e=="string"?t[e]||t[e.toLowerCase()]||e:t[e]||e}}g(Me,"type","icc"),g(Me,"multiSegment",!0),g(Me,"headerLength",18);const mi={4:J,8:function(i,e){return[i.getUint8(e),i.getUint8(e+1)>>4,i.getUint8(e+1)%16].map(t=>t.toString(10)).join(".")},12:J,16:J,20:J,24:function(i,e){const t=i.getUint16(e),n=i.getUint16(e+2)-1,r=i.getUint16(e+4),s=i.getUint16(e+6),a=i.getUint16(e+8),o=i.getUint16(e+10);return new Date(Date.UTC(t,n,r,s,a,o))},36:J,40:J,48:J,52:J,64:(i,e)=>i.getUint32(e),80:J};function J(i,e){return ce(i.getString(e,4))}R.set("icc",Me),I(D,"icc",[[4,"ProfileCMMType"],[8,"ProfileVersion"],[12,"ProfileClass"],[16,"ColorSpaceData"],[20,"ProfileConnectionSpace"],[24,"ProfileDateTime"],[36,"ProfileFileSignature"],[40,"PrimaryPlatform"],[44,"CMMFlags"],[48,"DeviceManufacturer"],[52,"DeviceModel"],[56,"DeviceAttributes"],[64,"RenderingIntent"],[68,"ConnectionSpaceIlluminant"],[80,"ProfileCreator"],[84,"ProfileID"],["Header","ProfileHeader"],["MS00","WCSProfiles"],["bTRC","BlueTRC"],["bXYZ","BlueMatrixColumn"],["bfd","UCRBG"],["bkpt","MediaBlackPoint"],["calt","CalibrationDateTime"],["chad","ChromaticAdaptation"],["chrm","Chromaticity"],["ciis","ColorimetricIntentImageState"],["clot","ColorantTableOut"],["clro","ColorantOrder"],["clrt","ColorantTable"],["cprt","ProfileCopyright"],["crdi","CRDInfo"],["desc","ProfileDescription"],["devs","DeviceSettings"],["dmdd","DeviceModelDesc"],["dmnd","DeviceMfgDesc"],["dscm","ProfileDescriptionML"],["fpce","FocalPlaneColorimetryEstimates"],["gTRC","GreenTRC"],["gXYZ","GreenMatrixColumn"],["gamt","Gamut"],["kTRC","GrayTRC"],["lumi","Luminance"],["meas","Measurement"],["meta","Metadata"],["mmod","MakeAndModel"],["ncl2","NamedColor2"],["ncol","NamedColor"],["ndin","NativeDisplayInfo"],["pre0","Preview0"],["pre1","Preview1"],["pre2","Preview2"],["ps2i","PS2RenderingIntent"],["ps2s","PostScript2CSA"],["psd0","PostScript2CRD0"],["psd1","PostScript2CRD1"],["psd2","PostScript2CRD2"],["psd3","PostScript2CRD3"],["pseq","ProfileSequenceDesc"],["psid","ProfileSequenceIdentifier"],["psvm","PS2CRDVMSize"],["rTRC","RedTRC"],["rXYZ","RedMatrixColumn"],["resp","OutputResponse"],["rhoc","ReflectionHardcopyOrigColorimetry"],["rig0","PerceptualRenderingIntentGamut"],["rig2","SaturationRenderingIntentGamut"],["rpoc","ReflectionPrintOutputColorimetry"],["sape","SceneAppearanceEstimates"],["scoe","SceneColorimetryEstimates"],["scrd","ScreeningDesc"],["scrn","Screening"],["targ","CharTarget"],["tech","Technology"],["vcgt","VideoCardGamma"],["view","ViewingConditions"],["vued","ViewingCondDesc"],["wtpt","MediaWhitePoint"]]);const Re={"4d2p":"Erdt Systems",AAMA:"Aamazing Technologies",ACER:"Acer",ACLT:"Acolyte Color Research",ACTI:"Actix Sytems",ADAR:"Adara Technology",ADBE:"Adobe",ADI:"ADI Systems",AGFA:"Agfa Graphics",ALMD:"Alps Electric",ALPS:"Alps Electric",ALWN:"Alwan Color Expertise",AMTI:"Amiable Technologies",AOC:"AOC International",APAG:"Apago",APPL:"Apple Computer",AST:"AST","AT&T":"AT&T",BAEL:"BARBIERI electronic",BRCO:"Barco NV",BRKP:"Breakpoint",BROT:"Brother",BULL:"Bull",BUS:"Bus Computer Systems","C-IT":"C-Itoh",CAMR:"Intel",CANO:"Canon",CARR:"Carroll Touch",CASI:"Casio",CBUS:"Colorbus PL",CEL:"Crossfield",CELx:"Crossfield",CGS:"CGS Publishing Technologies International",CHM:"Rochester Robotics",CIGL:"Colour Imaging Group, London",CITI:"Citizen",CL00:"Candela",CLIQ:"Color IQ",CMCO:"Chromaco",CMiX:"CHROMiX",COLO:"Colorgraphic Communications",COMP:"Compaq",COMp:"Compeq/Focus Technology",CONR:"Conrac Display Products",CORD:"Cordata Technologies",CPQ:"Compaq",CPRO:"ColorPro",CRN:"Cornerstone",CTX:"CTX International",CVIS:"ColorVision",CWC:"Fujitsu Laboratories",DARI:"Darius Technology",DATA:"Dataproducts",DCP:"Dry Creek Photo",DCRC:"Digital Contents Resource Center, Chung-Ang University",DELL:"Dell Computer",DIC:"Dainippon Ink and Chemicals",DICO:"Diconix",DIGI:"Digital","DL&C":"Digital Light & Color",DPLG:"Doppelganger",DS:"Dainippon Screen",DSOL:"DOOSOL",DUPN:"DuPont",EPSO:"Epson",ESKO:"Esko-Graphics",ETRI:"Electronics and Telecommunications Research Institute",EVER:"Everex Systems",EXAC:"ExactCODE",Eizo:"Eizo",FALC:"Falco Data Products",FF:"Fuji Photo Film",FFEI:"FujiFilm Electronic Imaging",FNRD:"Fnord Software",FORA:"Fora",FORE:"Forefront Technology",FP:"Fujitsu",FPA:"WayTech Development",FUJI:"Fujitsu",FX:"Fuji Xerox",GCC:"GCC Technologies",GGSL:"Global Graphics Software",GMB:"Gretagmacbeth",GMG:"GMG",GOLD:"GoldStar Technology",GOOG:"Google",GPRT:"Giantprint",GTMB:"Gretagmacbeth",GVC:"WayTech Development",GW2K:"Sony",HCI:"HCI",HDM:"Heidelberger Druckmaschinen",HERM:"Hermes",HITA:"Hitachi America",HP:"Hewlett-Packard",HTC:"Hitachi",HiTi:"HiTi Digital",IBM:"IBM",IDNT:"Scitex",IEC:"Hewlett-Packard",IIYA:"Iiyama North America",IKEG:"Ikegami Electronics",IMAG:"Image Systems",IMI:"Ingram Micro",INTC:"Intel",INTL:"N/A (INTL)",INTR:"Intra Electronics",IOCO:"Iocomm International Technology",IPS:"InfoPrint Solutions Company",IRIS:"Scitex",ISL:"Ichikawa Soft Laboratory",ITNL:"N/A (ITNL)",IVM:"IVM",IWAT:"Iwatsu Electric",Idnt:"Scitex",Inca:"Inca Digital Printers",Iris:"Scitex",JPEG:"Joint Photographic Experts Group",JSFT:"Jetsoft Development",JVC:"JVC Information Products",KART:"Scitex",KFC:"KFC Computek Components",KLH:"KLH Computers",KMHD:"Konica Minolta",KNCA:"Konica",KODA:"Kodak",KYOC:"Kyocera",Kart:"Scitex",LCAG:"Leica",LCCD:"Leeds Colour",LDAK:"Left Dakota",LEAD:"Leading Technology",LEXM:"Lexmark International",LINK:"Link Computer",LINO:"Linotronic",LITE:"Lite-On",Leaf:"Leaf",Lino:"Linotronic",MAGC:"Mag Computronic",MAGI:"MAG Innovision",MANN:"Mannesmann",MICN:"Micron Technology",MICR:"Microtek",MICV:"Microvitec",MINO:"Minolta",MITS:"Mitsubishi Electronics America",MITs:"Mitsuba",MNLT:"Minolta",MODG:"Modgraph",MONI:"Monitronix",MONS:"Monaco Systems",MORS:"Morse Technology",MOTI:"Motive Systems",MSFT:"Microsoft",MUTO:"MUTOH INDUSTRIES",Mits:"Mitsubishi Electric",NANA:"NANAO",NEC:"NEC",NEXP:"NexPress Solutions",NISS:"Nissei Sangyo America",NKON:"Nikon",NONE:"none",OCE:"Oce Technologies",OCEC:"OceColor",OKI:"Oki",OKID:"Okidata",OKIP:"Okidata",OLIV:"Olivetti",OLYM:"Olympus",ONYX:"Onyx Graphics",OPTI:"Optiquest",PACK:"Packard Bell",PANA:"Matsushita Electric Industrial",PANT:"Pantone",PBN:"Packard Bell",PFU:"PFU",PHIL:"Philips Consumer Electronics",PNTX:"HOYA",POne:"Phase One A/S",PREM:"Premier Computer Innovations",PRIN:"Princeton Graphic Systems",PRIP:"Princeton Publishing Labs",QLUX:"Hong Kong",QMS:"QMS",QPCD:"QPcard AB",QUAD:"QuadLaser",QUME:"Qume",RADI:"Radius",RDDx:"Integrated Color Solutions",RDG:"Roland DG",REDM:"REDMS Group",RELI:"Relisys",RGMS:"Rolf Gierling Multitools",RICO:"Ricoh",RNLD:"Edmund Ronald",ROYA:"Royal",RPC:"Ricoh Printing Systems",RTL:"Royal Information Electronics",SAMP:"Sampo",SAMS:"Samsung",SANT:"Jaime Santana Pomares",SCIT:"Scitex",SCRN:"Dainippon Screen",SDP:"Scitex",SEC:"Samsung",SEIK:"Seiko Instruments",SEIk:"Seikosha",SGUY:"ScanGuy.com",SHAR:"Sharp Laboratories",SICC:"International Color Consortium",SONY:"Sony",SPCL:"SpectraCal",STAR:"Star",STC:"Sampo Technology",Scit:"Scitex",Sdp:"Scitex",Sony:"Sony",TALO:"Talon Technology",TAND:"Tandy",TATU:"Tatung",TAXA:"TAXAN America",TDS:"Tokyo Denshi Sekei",TECO:"TECO Information Systems",TEGR:"Tegra",TEKT:"Tektronix",TI:"Texas Instruments",TMKR:"TypeMaker",TOSB:"Toshiba",TOSH:"Toshiba",TOTK:"TOTOKU ELECTRIC",TRIU:"Triumph",TSBT:"Toshiba",TTX:"TTX Computer Products",TVM:"TVM Professional Monitor",TW:"TW Casper",ULSX:"Ulead Systems",UNIS:"Unisys",UTZF:"Utz Fehlau & Sohn",VARI:"Varityper",VIEW:"Viewsonic",VISL:"Visual communication",VIVO:"Vivo Mobile Communication",WANG:"Wang",WLBR:"Wilbur Imaging",WTG2:"Ware To Go",WYSE:"WYSE Technology",XERX:"Xerox",XRIT:"X-Rite",ZRAN:"Zoran",Zebr:"Zebra Technologies",appl:"Apple Computer",bICC:"basICColor",berg:"bergdesign",ceyd:"Integrated Color Solutions",clsp:"MacDermid ColorSpan",ds:"Dainippon Screen",dupn:"DuPont",ffei:"FujiFilm Electronic Imaging",flux:"FluxData",iris:"Scitex",kart:"Scitex",lcms:"Little CMS",lino:"Linotronic",none:"none",ob4d:"Erdt Systems",obic:"Medigraph",quby:"Qubyx Sarl",scit:"Scitex",scrn:"Dainippon Screen",sdp:"Scitex",siwi:"SIWI GRAFIKA",yxym:"YxyMaster"},jt={scnr:"Scanner",mntr:"Monitor",prtr:"Printer",link:"Device Link",abst:"Abstract",spac:"Color Space Conversion Profile",nmcl:"Named Color",cenc:"ColorEncodingSpace profile",mid:"MultiplexIdentification profile",mlnk:"MultiplexLink profile",mvis:"MultiplexVisualization profile",nkpf:"Nikon Input Device Profile (NON-STANDARD!)"};I(j,"icc",[[4,Re],[12,jt],[40,Object.assign({},Re,jt)],[48,Re],[80,Re],[64,{0:"Perceptual",1:"Relative Colorimetric",2:"Saturation",3:"Absolute Colorimetric"}],["tech",{amd:"Active Matrix Display",crt:"Cathode Ray Tube Display",kpcd:"Photo CD",pmd:"Passive Matrix Display",dcam:"Digital Camera",dcpj:"Digital Cinema Projector",dmpc:"Digital Motion Picture Camera",dsub:"Dye Sublimation Printer",epho:"Electrophotographic Printer",esta:"Electrostatic Printer",flex:"Flexography",fprn:"Film Writer",fscn:"Film Scanner",grav:"Gravure",ijet:"Ink Jet Printer",imgs:"Photo Image Setter",mpfr:"Motion Picture Film Recorder",mpfs:"Motion Picture Film Scanner",offs:"Offset Lithography",pjtv:"Projection Television",rpho:"Photographic Paper Printer",rscn:"Reflective Scanner",silk:"Silkscreen",twax:"Thermal Wax Printer",vidc:"Video Camera",vidm:"Video Monitor"}]]);class Oe extends W{static canHandle(e,t,n){return e.getUint8(t+1)===237&&e.getString(t+4,9)==="Photoshop"&&this.containsIptc8bim(e,t,n)!==void 0}static headerLength(e,t,n){let r,s=this.containsIptc8bim(e,t,n);if(s!==void 0)return r=e.getUint8(t+s+7),r%2!=0&&(r+=1),r===0&&(r=4),s+8+r}static containsIptc8bim(e,t,n){for(let r=0;r<n;r++)if(this.isIptcSegmentHead(e,t+r))return r}static isIptcSegmentHead(e,t){return e.getUint8(t)===56&&e.getUint32(t)===943868237&&e.getUint16(t+4)===1028}parse(){let{raw:e}=this,t=this.chunk.byteLength-1,n=!1;for(let r=0;r<t;r++)if(this.chunk.getUint8(r)===28&&this.chunk.getUint8(r+1)===2){n=!0;let s=this.chunk.getUint16(r+3),a=this.chunk.getUint8(r+2),o=this.chunk.getLatin1String(r+5,s);e.set(a,this.pluralizeValue(e.get(a),o)),r+=4+s}else if(n)break;return this.translate(),this.output}pluralizeValue(e,t){return e!==void 0?e instanceof Array?(e.push(t),e):[e,t]:t}}g(Oe,"type","iptc"),g(Oe,"translateValues",!1),g(Oe,"reviveValues",!1),R.set("iptc",Oe),I(D,"iptc",[[0,"ApplicationRecordVersion"],[3,"ObjectTypeReference"],[4,"ObjectAttributeReference"],[5,"ObjectName"],[7,"EditStatus"],[8,"EditorialUpdate"],[10,"Urgency"],[12,"SubjectReference"],[15,"Category"],[20,"SupplementalCategories"],[22,"FixtureIdentifier"],[25,"Keywords"],[26,"ContentLocationCode"],[27,"ContentLocationName"],[30,"ReleaseDate"],[35,"ReleaseTime"],[37,"ExpirationDate"],[38,"ExpirationTime"],[40,"SpecialInstructions"],[42,"ActionAdvised"],[45,"ReferenceService"],[47,"ReferenceDate"],[50,"ReferenceNumber"],[55,"DateCreated"],[60,"TimeCreated"],[62,"DigitalCreationDate"],[63,"DigitalCreationTime"],[65,"OriginatingProgram"],[70,"ProgramVersion"],[75,"ObjectCycle"],[80,"Byline"],[85,"BylineTitle"],[90,"City"],[92,"Sublocation"],[95,"State"],[100,"CountryCode"],[101,"Country"],[103,"OriginalTransmissionReference"],[105,"Headline"],[110,"Credit"],[115,"Source"],[116,"CopyrightNotice"],[118,"Contact"],[120,"Caption"],[121,"LocalCaption"],[122,"Writer"],[125,"RasterizedCaption"],[130,"ImageType"],[131,"ImageOrientation"],[135,"LanguageIdentifier"],[150,"AudioType"],[151,"AudioSamplingRate"],[152,"AudioSamplingResolution"],[153,"AudioDuration"],[154,"AudioOutcue"],[184,"JobID"],[185,"MasterDocumentID"],[186,"ShortDocumentID"],[187,"UniqueDocumentID"],[188,"OwnerID"],[200,"ObjectPreviewFileFormat"],[201,"ObjectPreviewFileVersion"],[202,"ObjectPreviewData"],[221,"Prefs"],[225,"ClassifyState"],[228,"SimilarityIndex"],[230,"DocumentNotes"],[231,"DocumentHistory"],[232,"ExifCameraInfo"],[255,"CatalogSets"]]),I(j,"iptc",[[10,{0:"0 (reserved)",1:"1 (most urgent)",2:"2",3:"3",4:"4",5:"5 (normal urgency)",6:"6",7:"7",8:"8 (least urgent)",9:"9 (user-defined priority)"}],[75,{a:"Morning",b:"Both Morning and Evening",p:"Evening"}],[131,{L:"Landscape",P:"Portrait",S:"Square"}]]);let ot=null;async function _i(){if(!ot)try{const n=(await import("./joraw-1Lq5hXK7.js")).default;if(typeof n!="function")throw new Error("JoRaw WASM import failed");const r=new URL("/assets/joraw-DraTMNgX.wasm",import.meta.url).href;ot=n({locateFile:(s,a)=>s.endsWith("joraw.wasm")?r:a+s})}catch(t){throw console.error("Failed to load joraw.js:",t),t}const i=await ot,e=i.LibRaw||i.JoRaw;if(!e)throw new Error("JoRaw class not found");return e}const yi=async i=>{var n,r,s,a,o,d,u,c,l,p,h,m,_,b;const e=await _i(),t=new e;try{if(await t.open(i,{}),typeof t.getRawImage!="function")throw new Error("WASM mismatch");const w=t.getRawImage();let C=new Uint16Array(w.data);const y=await t.metadata(!0);let k={...y};try{const f=await pi.parse(i.buffer);f&&(k={...k,...f})}catch(f){console.warn("exifr parsing failed for RAW buffer",f)}const S=((n=y.idata)==null?void 0:n.filters)||0,A=((r=y.idata)==null?void 0:r.colors)||0,F=S===0&&A===3,X=S===9;let x=[0,0,0,0],T=!1,V=null;if(t.getBlackLevels)try{const f=t.getBlackLevels();V=f,f.dng_cblack&&f.dng_cblack.length===4&&Array.from(f.dng_cblack).some(L=>L>0)?(x=Array.from(f.dng_cblack).map(Number),T=!0):f.cblack&&f.cblack.length===4&&Array.from(f.cblack).some(L=>L>0)?(x=Array.from(f.cblack).map(Number),T=!0):typeof f.black=="number"&&f.black>0&&(x=[f.black,f.black,f.black,f.black],T=!0)}catch(f){console.warn("getBlackLevels binding failed",f)}if(!T){let f=[];if((s=y.color_data)!=null&&s.cblack_rawpy_style)f=y.color_data.cblack_rawpy_style;else if((o=(a=y.color_data)==null?void 0:a.dng_levels)!=null&&o.dng_cblack)f=y.color_data.dng_levels.dng_cblack;else if(((d=y.black_level_per_channel)==null?void 0:d.length)>=4)f=y.black_level_per_channel;else if(((u=y.cblack)==null?void 0:u.length)>=4)f=y.cblack;else if(((l=(c=y.color)==null?void 0:c.cblack)==null?void 0:l.length)>=4)f=y.color.cblack;else{const L=y.black_level||y.color_data&&y.color_data.black||0;f=[L,L,L,L]}x=[Number(f[0])||0,Number(f[1])||0,Number(f[2])||0,Number(f[3])||0]}const z=(...f)=>{for(const L of f){if(Array.isArray(L)||ArrayBuffer.isView(L)){const M=Math.max(...Array.from(L).map(Number).filter(P=>Number.isFinite(P)&&P>0));if(Number.isFinite(M)&&M>0)return M;continue}const U=Number(L);if(Number.isFinite(U)&&U>0)return U}return null},K=Number.isFinite(Number(w.bits))&&Number(w.bits)>0?Math.pow(2,Number(w.bits))-1:null,G=z(y.white_level,(h=(p=y.color_data)==null?void 0:p.dng_levels)==null?void 0:h.dng_whitelevel,(m=y.color_data)==null?void 0:m.maximum,V==null?void 0:V.maximum,(_=y.color_data)==null?void 0:_.fmaximum,(b=y.color_data)==null?void 0:b.data_maximum,K)||16383;return{data:C,width:w.width,height:w.height,bayerPattern:y.color_desc||"RGGB",blackLevels:x,whiteLevel:G,metadata:k,isThreePlane:F,threePlaneTransfer:F?"linear":void 0,isXTrans:X}}finally{t.delete?t.delete():t.close()}};async function wi(i){if(Rn(i)){const n=await Nn(i);if(!n)throw new Error("Sony cRAW HQ decoder did not return image data.");return n.rawImageData}if(Wn(i)){const n=await jn(i);if(!n)throw new Error("Nikon HE decoder did not return image data.");return n}const t=new Uint8Array(i);return yi(t)}function bi(i,e,t,n,r,s=!1,a=!1,o={rows:!1,rowGroup:0,cols:!1,colGroup:0},d=!0,u=1,c=!1){const l=Math.max(0,Math.floor(n.x)),p=Math.max(0,Math.floor(n.y)),h=Math.min(e,Math.ceil(n.x+n.w)),m=Math.min(t,Math.ceil(n.y+n.h)),_=h-l,b=m-p;if(_<=1||b<=1)return null;const w=c&&i.length>=e*t*3;if(w){if(r<0||r>2)return null}else if(r<0||r>3)return null;let C=0,y=0;if(w)C=h-l,y=m-p;else{for(let f=l;f<h;f++)f%2===(r&1)&&C++;for(let f=p;f<m;f++)f%2===r>>1&&y++}if(C<4||y<4)return null;const k=Math.max(1,u),S=Math.floor(2*C/(k+1)),A=Math.floor(2*y/(k+1));if(S<4||A<4)return null;const F=new Fe(S),X=new Fe(A),x=new Float32Array(Math.floor(S/2)+1),T=new Float32Array(Math.floor(A/2)+1);let V=0,z=0;const K=new Float32Array(C),G=new Float32Array(y);for(let f=p;f<m;f++){const L=f%2;if(!w&&r>>1!==L||o.rows&&(w?f&1:f>>1&1)!==o.rowGroup)continue;const U=f*e;let M=0;for(let P=l;P<h;P++){if(!w&&P%2!==(r&1)||M>=C)continue;const q=w?i[(U+P)*3+r]:i[U+P],qe=Number.isFinite(q)?q:0;K[M++]=qe}for(let P=0;P<k;P++){const q=Math.floor(P*S/2);if(q+S>C)break;F.calculateSpectrumWindow(K,q,x,a),V++}if(!d)break}for(let f=l;f<h;f++){const L=f%2;if(!w&&(r&1)!==L||o.cols&&(w?f&1:f>>1&1)!==o.colGroup)continue;let U=0;for(let M=p;M<m;M++){if(!w&&M%2!==r>>1||U>=y)continue;const P=w?i[(M*e+f)*3+r]:i[M*e+f],q=Number.isFinite(P)?P:0;G[U++]=q}for(let M=0;M<k;M++){const P=Math.floor(M*A/2);if(P+A>y)break;X.calculateSpectrumWindow(G,P,T,a),z++}if(!d)break}if(V>0)for(let f=0;f<x.length;f++)x[f]/=V;if(z>0)for(let f=0;f<T.length;f++)T[f]/=z;return{avgHoriz:Array.from(x),avgVert:Array.from(T),sampleCount:{h:V,v:z}}}const Xt=192*1024*1024,ze=new Map;let De=0,Kt=Promise.resolve();const xi=i=>{i.lastAccess=Date.now()},ki=(i=new Set)=>{if(De<=Xt)return;const e=Array.from(ze.entries()).filter(([t])=>!i.has(t)).sort((t,n)=>t[1].lastAccess-n[1].lastAccess);for(const[t,n]of e){if(De<=Xt)break;De-=n.sizeBytes,ze.delete(t)}},Si=i=>{var e;return i.data.byteLength+(((e=i.floatData)==null?void 0:e.byteLength)||0)},un=async i=>{const e=ze.get(i.key);if(e)return xi(e),e.raw;if(!i.file)throw new Error(`Worker cache miss for ${i.key}`);const t=await i.file.arrayBuffer(),n=await wi(t),r={raw:n,sizeBytes:Si(n),lastAccess:Date.now()};return ze.set(i.key,r),De+=r.sizeBytes,ki(new Set([i.key])),n},vi=async i=>{for(const e of i.files)await un(e);self.postMessage({id:i.id,type:"result",success:!0})},Ci=async i=>{const e=new Map;for(const r of i.files){const s=await un(r);for(const a of i.channelIndices){const o=bi(s.data,s.width,s.height,i.roi,a,!1,i.useHanning,i.interlace,i.doAverage,i.nSegments,!!s.isThreePlane);if(!o)continue;let d=e.get(a);if(!d)d={sumH:new Float32Array(o.avgHoriz.length),sumV:new Float32Array(o.avgVert.length),sampleCount:o.sampleCount,contributingFrames:0},e.set(a,d);else if(d.sumH.length!==o.avgHoriz.length||d.sumV.length!==o.avgVert.length)throw new Error("Stack FFT requires matching frame dimensions and channel layout.");for(let u=0;u<o.avgHoriz.length;u++)d.sumH[u]+=o.avgHoriz[u];for(let u=0;u<o.avgVert.length;u++)d.sumV[u]+=o.avgVert[u];d.sampleCount=o.sampleCount,d.contributingFrames+=1}self.postMessage({id:i.id,type:"progress",completed:1})}const t=Array.from(e.entries()).map(([r,s])=>({ch:r,sumH:s.sumH,sumV:s.sumV,sampleCount:s.sampleCount,contributingFrames:s.contributingFrames})),n=[];for(const r of t)n.push(r.sumH.buffer,r.sumV.buffer);self.postMessage({id:i.id,type:"result",success:!0,processedFrames:i.files.length,results:t},n)};self.onmessage=i=>{const e=i.data;Kt=Kt.then(async()=>{try{e.type==="preload"?await vi(e):await Ci(e)}catch(t){self.postMessage({id:e.id,type:"result",success:!1,error:(t==null?void 0:t.message)||"Power spectrum worker failed"})}})};
