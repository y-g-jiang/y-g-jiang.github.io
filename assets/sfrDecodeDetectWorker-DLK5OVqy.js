var Ns=Object.defineProperty;var Rs=(n,t,e)=>t in n?Ns(n,t,{enumerable:!0,configurable:!0,writable:!0,value:e}):n[t]=e;var bt=(n,t,e)=>Rs(n,typeof t!="symbol"?t+"":t,e);function Es(n){const[[t,e,i],[r,s,o],[a,l,u]]=n,d=t*(s*u-o*l)-e*(r*u-o*a)+i*(r*l-s*a);if(Math.abs(d)<1e-10)return null;const c=1/d;return[[(s*u-o*l)*c,(i*l-e*u)*c,(e*o-i*s)*c],[(o*a-r*u)*c,(t*u-i*a)*c,(i*r-t*o)*c],[(r*l-s*a)*c,(e*a-t*l)*c,(t*s-e*r)*c]]}function Ls(n,t){const[[e,i,r],[s,o,a],[l,u,d]]=n,[c,g,h]=t;return[e*c+i*g+r*h,s*c+o*g+a*h,l*c+u*g+d*h]}class Cn{constructor(t){bt(this,"size");bt(this,"isPowerOfTwo");bt(this,"_real");bt(this,"_imag");bt(this,"_scratch",null);bt(this,"_rev",null);bt(this,"_m",0);bt(this,"_internalFFT",null);bt(this,"_chirpReal",null);bt(this,"_chirpImag",null);bt(this,"_bReal",null);bt(this,"_bImag",null);bt(this,"_hanning",null);bt(this,"_windowSumSq",0);this.size=t,this.isPowerOfTwo=(t&t-1)===0&&t>0,this._real=new Float32Array(t),this._imag=new Float32Array(t),this.isPowerOfTwo?this.initRadix2():this.initBluestein()}initRadix2(){const t=this.size,e=Math.log2(t);this._rev=new Uint32Array(t);for(let i=0;i<t;i++){let r=0,s=i;for(let o=0;o<e;o++)r=r<<1|s&1,s>>>=1;this._rev[i]=r}}initBluestein(){const t=this.size;this._m=Math.pow(2,Math.ceil(Math.log2(2*t-1))),this._internalFFT=new Cn(this._m),this._chirpReal=new Float32Array(t),this._chirpImag=new Float32Array(t);for(let r=0;r<t;r++){const s=-Math.PI*(r*r)/t;this._chirpReal[r]=Math.cos(s),this._chirpImag[r]=Math.sin(s)}const e=new Float32Array(this._m),i=new Float32Array(this._m);for(let r=0;r<t;r++)e[r]=this._chirpReal[r],i[r]=-this._chirpImag[r];for(let r=1;r<t;r++)e[this._m-r]=e[r],i[this._m-r]=i[r];this._internalFFT.transform(e,i),this._bReal=new Float32Array(this._internalFFT._real),this._bImag=new Float32Array(this._internalFFT._imag)}initHanning(){if(this._hanning)return;const t=this.size;this._hanning=new Float32Array(t);let e=0;for(let i=0;i<t;i++){const r=.5*(1-Math.cos(2*Math.PI*i/(t-1)));this._hanning[i]=r,e+=r*r}this._windowSumSq=e}transform(t,e){this.isPowerOfTwo?this.transformRadix2(t,e):this.transformBluestein(t,e)}transformRadix2(t,e){const i=this.size,r=this._rev,s=this._real,o=this._imag;if(t===s)for(let a=0;a<i;a++){const l=r[a];if(a<l){const u=s[a],d=o[a];s[a]=s[l],o[a]=o[l],s[l]=u,o[l]=d}}else for(let a=0;a<i;a++){const l=r[a];s[a]=t[l],o[a]=e?e[l]:0}for(let a=2;a<=i;a*=2){const l=a/2,u=-2*Math.PI/a,d=Math.cos(u),c=Math.sin(u);for(let g=0;g<i;g+=a){let h=1,f=0;for(let p=0;p<l;p++){const m=g+p,y=g+p+l,x=h*s[y]-f*o[y],b=h*o[y]+f*s[y],M=s[m],C=o[m];s[m]=M+x,o[m]=C+b,s[y]=M-x,o[y]=C-b;const S=h*d-f*c,w=h*c+f*d;h=S,f=w}}}}transformBluestein(t,e){const i=this.size,r=this._m,s=this._internalFFT,o=s._real,a=s._imag;o.fill(0),a.fill(0);for(let c=0;c<i;c++){const g=t[c],h=e?e[c]:0,f=this._chirpReal[c],p=this._chirpImag[c];o[c]=g*f-h*p,a[c]=g*p+h*f}s.transformRadix2(o,a);for(let c=0;c<r;c++){const g=s._real[c],h=s._imag[c],f=this._bReal[c],p=this._bImag[c];s._real[c]=g*f-h*p,s._imag[c]=g*p+h*f}const l=s._real,u=s._imag;for(let c=0;c<r;c++)u[c]=-u[c];s.transformRadix2(l,u);const d=1/r;for(let c=0;c<i;c++){const g=s._real[c]*d,h=-s._imag[c]*d,f=this._chirpReal[c],p=this._chirpImag[c];this._real[c]=g*f-h*p,this._imag[c]=g*p+h*f}}calculateSpectrum(t,e,i=!1){const r=this.size;let s=0;for(let d=0;d<r;d++)s+=t[d];const o=s/r;this._scratch||(this._scratch=new Float32Array(r));const a=this._scratch;if(i){this.initHanning();const d=this._hanning;for(let c=0;c<r;c++)a[c]=(t[c]-o)*d[c]}else for(let d=0;d<r;d++)a[d]=t[d]-o;this.transform(a);const l=e.length;let u=1/r;i&&this._windowSumSq>0&&(u=1/this._windowSumSq);for(let d=0;d<l;d++){const c=this._real[d],g=this._imag[d];e[d]+=(c*c+g*g)*u}}calculateSpectrumWindow(t,e,i,r=!1){const s=this.size;let o=0;for(let c=0;c<s;c++)o+=t[e+c];const a=o/s;this._scratch||(this._scratch=new Float32Array(s));const l=this._scratch;if(r){this.initHanning();const c=this._hanning;for(let g=0;g<s;g++)l[g]=(t[e+g]-a)*c[g]}else for(let c=0;c<s;c++)l[c]=t[e+c]-a;this.transform(l);const u=i.length;let d=1/s;r&&this._windowSumSq>0&&(d=1/this._windowSumSq);for(let c=0;c<u;c++){const g=this._real[c],h=this._imag[c];i[c]+=(g*g+h*h)*d}}}const Ds={"Sony ILCE-7RM5":"0.82 -0.2976 -0.0719 -0.4296 1.2053 0.2532 -0.0429 0.1282 0.5774"};let vi=null;async function Bs(n){return vi||(vi=(async()=>{if(typeof window.loadPyodide!="function")throw new Error("Pyodide missing: window.loadPyodide not found.");const t=await window.loadPyodide();return await t.loadPackage("numpy"),t})()),vi}var Us=`#!/usr/bin/env python3
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
`,Os=`#!/usr/bin/env python3
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

from llvc3_bitstream_probe import RAW_STREAM_OFFSET, find_raw_subifd


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


def load_packet(arw: Path, group: int, index: int) -> tuple[bytes, dict]:
    from llvc3_bitstream_probe import parse_directory, parse_llvc_header, parse_packet

    raw_info, strip = find_raw_subifd(arw)
    stream = strip[RAW_STREAM_OFFSET:]
    parse_llvc_header(stream)
    _consumed, entries, _groups = parse_directory(stream)
    entry = next(e for e in entries if e.group == group and e.index == index)
    packet_info = parse_packet(stream, entry)
    packet = stream[entry.start : entry.start + entry.length]
    info = json.loads(json.dumps(packet_info, default=lambda o: o.__dict__))
    info["raw_width"] = raw_info.width
    info["raw_height"] = raw_info.height
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


def decode_packet_components(arw: Path, group: int, index: int, out_prefix: Path | None = None) -> dict[str, object]:
    """Decode a type-1 or type-3 packet into one or three int32 component arrays."""

    packet, info = load_packet(arw, group, index)
    packet_type = info["type2"]
    components = 1 if packet_type == 1 else 3
    width = infer_packet_width(group, packet_type, int(info.get("raw_width", 7040)))
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


def decode_packet_arrays(arw: Path, group: int, index: int) -> tuple[list["object"], dict[str, object]]:
    """Decode a packet and return its int32 component arrays in memory."""

    import numpy as np

    packet, info = load_packet(arw, group, index)
    packet_type = info["type2"]
    components = 1 if packet_type == 1 else 3
    width = infer_packet_width(group, packet_type, int(info.get("raw_width", 7040)))
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
        "shape": [int(arrays[0].shape[0]), int(arrays[0].shape[1])],
        "row_multiplier": row_multiplier,
        "components": components,
        "final_state": final_states[-1] if final_states else [],
    }
    return arrays, meta


def decode_type1_packet(arw: Path, group: int, index: int, out: Path | None = None) -> dict[str, object]:
    """Decode an independently parsed type-1 packet into int32 coefficient rows."""

    packet, info = load_packet(arw, group, index)
    if info["type2"] != 1:
        raise ValueError(f"packet g{group}i{index} is type {info['type2']}, not type 1")
    width = infer_packet_width(group, info["type2"], int(info.get("raw_width", 7040)))
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
`,Xs=`#!/usr/bin/env python3
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

    def edge_fn(x: np.ndarray) -> np.ndarray:
        xi = x.astype(np.int32)
        signed_half_step = np.where(xi > 0, 1, np.where(xi < 0, -1, 0))
        if edge_mode == "even":
            mask = ((xi & 1) == 0) & (xi != 0)
        elif edge_mode == "odd":
            mask = (xi & 1) != 0
        else:
            raise ValueError(f"unknown LLVC3 edge_mode {edge_mode!r}")
        return (2 * xi + np.where(mask, signed_half_step, 0)).astype(np.int32)

    hh[:edge_rows] = edge_fn(sub1[:edge_rows])
    hh[edge_rows : h - edge_rows] = sub2[np.arange(edge_rows, h - edge_rows) + edge_rows]
    bottom_hh = edge_fn(sub1[h : h + edge_rows])
    if bottom_hh_extra is not None:
        extra = np.asarray(bottom_hh_extra, dtype=np.int32)
        if extra.shape != bottom_hh.shape:
            raise ValueError(f"bottom_hh_extra shape {extra.shape} != bottom edge {bottom_hh.shape}")
        bottom_hh = bottom_hh + extra
    hh[h - edge_rows :] = bottom_hh

    low_horizontal = sony_inv53_1d(ll_i, lh, axis=0)
    high_horizontal = sony_inv53_1d(hl, hh, axis=0)
    return sony_inv53_1d(low_horizontal, high_horizontal, axis=1)


def synthesize_llvc3_final_green(ll: np.ndarray, detail: np.ndarray) -> np.ndarray:
    """Final CFA-green reconstruction, from the 0x1ab570 path.

    Not the same 2-D 5/3 inverse used by groups 1..3. It expands the half-width
    green lowpass into both RGGB green sites. The row offsets look odd because
    the native line buffer keeps four guard rows at the top and a few latency
    rows at the bottom.
    """

    ll_i = np.asarray(ll, dtype=np.int32)
    det = np.asarray(detail, dtype=np.int32)
    h, w = ll_i.shape
    if det.shape[1] != w or det.shape[0] < h + 4:
        raise ValueError(f"unexpected final green shapes: ll={ll_i.shape}, detail={det.shape}")

    selected = np.empty((h, w), dtype=np.int32)
    top = min(4, h)
    selected[:top] = det[:top]
    if h > 4:
        selected[4:] = det[8 : 8 + (h - 4)]

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
        v0_color = average(final_green_pair) + 2 * residual
    """

    g = np.asarray(v1_green, dtype=np.int32)
    r = np.asarray(v1_red, dtype=np.int32)
    b = np.asarray(v1_blue, dtype=np.int32)
    fg = np.asarray(full_green, dtype=np.int32)
    if fg.shape != (g.shape[0], g.shape[1] * 2):
        raise ValueError(f"unexpected final green shape: v1={g.shape}, full={fg.shape}")
    if r.shape != g.shape or b.shape != g.shape:
        raise ValueError(f"unexpected v1 color shapes: green={g.shape}, red={r.shape}, blue={b.shape}")

    gavg = (fg[:, 0::2] + fg[:, 1::2]) >> 1
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
`;const dr=512,Gs=Ds["Sony ILCE-7RM5"].split(/\s+/).map(Number).filter(Number.isFinite);let Pi=null;function Ys(n){if(n.byteLength<8)return null;const t=n.getUint16(0,!1);return t===18761?!0:t===19789?!1:null}function zr(n,t,e){const i=Math.min(n.length,t+e);let r="";for(let s=t;s<i;s++){const o=n[s];if(o===0)break;r+=String.fromCharCode(o)}return r.trim()}function Ht(n,t,e,i,r){const s=e===1||e===2||e===7?1:e===3||e===8?2:e===4||e===9?4:0;if(!s)return[];const o=s*i,a=o<=4?r:n.getUint32(r,t);if(a<0||a+o>n.byteLength)return[];const l=[];for(let u=0;u<i;u++){const d=a+u*s;e===1||e===2||e===7?l.push(n.getUint8(d)):e===3?l.push(n.getUint16(d,t)):e===8?l.push(n.getInt16(d,t)):e===4?l.push(n.getUint32(d,t)):e===9&&l.push(n.getInt32(d,t))}return l}function fr(n,t,e,i,r,s){if(i!==2||r<=0)return"";const o=r<=4?s:t.getUint32(s,e);return o<0||o>=n.length?"":zr(n,o,r)}function Vr(n){const t=new Uint8Array(n),e=new DataView(n),i=Ys(e);if(i===null||e.getUint16(2,i)!==42)return null;const s=c=>e.getUint16(c,i),o=c=>e.getUint32(c,i),a=[o(4)],l=new Set;let u="",d="";for(;a.length;){const c=a.pop();if(l.has(c)||c<=0||c+2>e.byteLength)continue;l.add(c);const g=s(c);if(c+2+g*12+4>e.byteLength)continue;const h=new Map;for(let M=0;M<g;M++){const C=c+2+M*12,S=s(C),w=s(C+2),P=o(C+4),T=C+8;h.set(S,{type:w,count:P,valueOffset:T})}const f=h.get(271),p=h.get(272);f&&!u&&(u=fr(t,e,i,f.type,f.count,f.valueOffset)),p&&!d&&(d=fr(t,e,i,p.type,p.count,p.valueOffset));const m=h.get(330);if(m){const M=Ht(e,i,m.type,m.count,m.valueOffset);for(const C of M)a.push(C)}const y=h.get(259),x=h.get(262);if(y&&x){const M=Ht(e,i,y.type,y.count,y.valueOffset)[0],C=Ht(e,i,x.type,x.count,x.valueOffset)[0];if(M===32766&&C===32803){const S=Ht(e,i,h.get(256).type,h.get(256).count,h.get(256).valueOffset)[0],w=Ht(e,i,h.get(257).type,h.get(257).count,h.get(257).valueOffset)[0],P=Ht(e,i,h.get(258).type,h.get(258).count,h.get(258).valueOffset)[0],T=Ht(e,i,h.get(273).type,h.get(273).count,h.get(273).valueOffset)[0],F=Ht(e,i,h.get(279).type,h.get(279).count,h.get(279).valueOffset)[0],k=h.get(33422)?Ht(e,i,h.get(33422).type,h.get(33422).count,h.get(33422).valueOffset):[0,1,1,2],_=h.get(29456)?Ht(e,i,h.get(29456).type,h.get(29456).count,h.get(29456).valueOffset):[],v=h.get(50717)?Ht(e,i,h.get(50717).type,h.get(50717).count,h.get(50717).valueOffset)[0]:16383,I=h.get(50719)?Ht(e,i,h.get(50719).type,h.get(50719).count,h.get(50719).valueOffset):[],E=h.get(50720)?Ht(e,i,h.get(50720).type,h.get(50720).count,h.get(50720).valueOffset):[];if(T+dr+16>t.length||T+F>t.length)return null;const L=T+dr,z=zr(t,L,4),U=t[L+8]<<8|t[L+9],R=t[L+10]<<8|t[L+11],O=t[L+12]<<8|t[L+13],B=t[L+14]<<8|t[L+15],Y=O>>4&63,Z=B>>13,X=B>>10&3;if(z!=="A000"||U!==S||R*2!==w||Y!==16||Z!==3||X!==3)return null;const J=[Number(_[0]??512),Number(_[1]??_[0]??512),Number(_[2]??_[0]??512),Number(_[3]??_[0]??512)],W=k.slice(0,4).map(tt=>tt===0?"R":tt===2?"B":"G").join("")||"RGGB";return{width:S,height:w,bitsPerSample:P,compression:M,photometric:C,blackLevel:J,whiteLevel:Number(v||16383),cfaPattern:W,defaultCropOrigin:I.length>=2?[Number(I[0]),Number(I[1])]:void 0,defaultCropSize:E.length>=2?[Number(E[0]),Number(E[1])]:void 0,make:u||"SONY",model:d||"ILCE-7M5"}}}const b=o(c+2+g*12);b&&a.push(b)}return null}async function zs(n){return Pi||(Pi=(async()=>{const t=await Bs();return t.__jtrSonyCrawHqDecoderReady||(await t.FS.mkdirTree("/sony_craw_hq"),await t.FS.writeFile("/sony_craw_hq/llvc3_bitstream_probe.py",Us),await t.FS.writeFile("/sony_craw_hq/llvc3_entropy.py",Os),await t.FS.writeFile("/sony_craw_hq/llvc3_math.py",Xs),await t.runPythonAsync(`
import sys
if "/sony_craw_hq" not in sys.path:
    sys.path.insert(0, "/sony_craw_hq")
from pathlib import Path
import numpy as np
from llvc3_bitstream_probe import find_raw_subifd, parse_llvc_header, RAW_STREAM_OFFSET
from llvc3_entropy import decode_packet_arrays, integrate_type1_coefficients
from llvc3_math import apply_sample_lut, clamp_signed_to_code_range, recombine_rggb, signed_to_sample
from llvc3_math import finalize_llvc3_color_planes, synthesize_llvc3_final_green, synthesize_llvc3_level_stride

def jtr_decode_signed_planes(arw):
    raw_info, _strip = find_raw_subifd(arw)
    low_rows = raw_info.height // 16

    g0, _meta = decode_packet_arrays(arw, 0, 0)
    green = integrate_type1_coefficients(g0[0][:low_rows], 2048) - 2048

    r0, _meta = decode_packet_arrays(arw, 0, 1)
    red_residual = integrate_type1_coefficients(r0[0][:low_rows], 0)

    b0, _meta = decode_packet_arrays(arw, 0, 2)
    blue_residual = integrate_type1_coefficients(b0[0][:low_rows], 0)

    for group, edge_rows in ((1, 0), (2, 1), (3, 2)):
        old_green = green
        old_red_residual = red_residual
        old_blue_residual = blue_residual

        planes, _meta = decode_packet_arrays(arw, group, 0)
        green = synthesize_llvc3_level_stride(old_green, planes[0], planes[1], planes[2], edge_rows)

        planes, _meta = decode_packet_arrays(arw, group, 1)
        edge_mode = "odd" if group == 3 else "even"
        red_residual = synthesize_llvc3_level_stride(
            old_red_residual, planes[0], planes[1], planes[2], edge_rows, edge_mode=edge_mode
        )

        planes, _meta = decode_packet_arrays(arw, group, 2)
        blue_residual = synthesize_llvc3_level_stride(
            old_blue_residual, planes[0], planes[1], planes[2], edge_rows, edge_mode=edge_mode
        )

    g4, _meta = decode_packet_arrays(arw, 4, 0)
    full_green = synthesize_llvc3_final_green(green, g4[0])
    v1_red = green + 2 * red_residual
    v1_blue = green + 2 * blue_residual
    return finalize_llvc3_color_planes(green, v1_red, v1_blue, full_green)

def jtr_decode_sony_craw_hq(arw_bytes, lut_bytes=None):
    path = Path("/tmp/jtr_sony_craw_hq_input.arw")
    path.write_bytes(bytes(arw_bytes))
    raw_info, strip = find_raw_subifd(path)
    header = parse_llvc_header(strip[RAW_STREAM_OFFSET:])
    if header.coded_width != raw_info.width or header.logical_height != raw_info.height or header.component_count != 3:
        raise ValueError(f"unexpected ARW6/LLVC3 header/raw mismatch: raw={raw_info}, header={header}")
    if raw_info.width % 16 or raw_info.height % 16:
        raise ValueError(f"decoder expects dimensions divisible by 16, got {raw_info.width}x{raw_info.height}")
    signed_c0, signed_c1, signed_c2 = jtr_decode_signed_planes(path)
    if lut_bytes is not None:
        lut = np.frombuffer(bytes(lut_bytes), dtype="<u2")
        if lut.size:
            if lut.size < 65536:
                lut = np.pad(lut, (0, 65536 - lut.size), constant_values=int(lut[-1]))
            lut = lut[:65536].astype(np.uint16)
            sample_c0 = apply_sample_lut(signed_to_sample(clamp_signed_to_code_range(signed_c0)), lut)
            sample_c1 = apply_sample_lut(signed_to_sample(clamp_signed_to_code_range(signed_c1)), lut)
            sample_c2 = apply_sample_lut(signed_to_sample(clamp_signed_to_code_range(signed_c2)), lut)
        else:
            sample_c0 = signed_to_sample(signed_c0)
            sample_c1 = signed_to_sample(signed_c1)
            sample_c2 = signed_to_sample(signed_c2)
    else:
        sample_c0 = signed_to_sample(signed_c0)
        sample_c1 = signed_to_sample(signed_c1)
        sample_c2 = signed_to_sample(signed_c2)
    raw = recombine_rggb(sample_c0, sample_c1, sample_c2)
    return raw.astype("<u2", copy=False).tobytes()
`),t.__jtrSonyCrawHqDecoderReady=!0),t})()),Pi}function Vs(n){return Vr(n)}async function js(n,t){const e=typeof performance<"u"?performance.now():Date.now(),i=Vr(n);if(!i)return null;const r=await zs(),s=typeof performance<"u"?performance.now():Date.now(),o=new Uint8Array(n),a=await fetch(new URL("/assets/sony_llvc3_static_lut4096_padded_u16-FsVBk-IV.bin",import.meta.url));if(!a.ok)throw new Error(`Failed to load Sony LLVC3 sample LUT: HTTP ${a.status}`);const l=new Uint8Array(await a.arrayBuffer()),u=typeof performance<"u"?performance.now():Date.now();r.globals.set("jtr_sony_arw_bytes",o),r.globals.set("jtr_sony_lut_bytes",l);const d=await r.runPythonAsync("jtr_decode_sony_craw_hq(jtr_sony_arw_bytes.to_py(), jtr_sony_lut_bytes.to_py())"),c=typeof performance<"u"?performance.now():Date.now(),g=d.toJs();typeof d.destroy=="function"&&d.destroy(),r.globals.delete("jtr_sony_arw_bytes"),r.globals.delete("jtr_sony_lut_bytes");const h=new Uint8Array(g.byteLength);h.set(g);const f=new Uint16Array(h.buffer),p=typeof performance<"u"?performance.now():Date.now();if(f.length!==i.width*i.height)throw new Error(`Sony cRAW HQ decoded size mismatch: got ${f.length}, expected ${i.width*i.height}`);const m=i.model||"ILCE-7M5",y=m.startsWith("Sony ")?m:`Sony ${m}`,x=m==="ILCE-7M5"?Gs:null,b={make:i.make||"SONY",model:m,camera_make:i.make||"SONY",camera_model:m,UniqueCameraModel:y,sourceFormat:"Sony cRAW HQ / LLVC3",sonyCrawHq:i,color_desc:i.cfaPattern,black_level_per_channel:i.blackLevel,white_level:i.whiteLevel,color_matrix:x&&x.length===9?x:void 0,idata:{filters:2492765332,colors:3},color_data:{cblack_rawpy_style:i.blackLevel,dng_levels:{dng_cblack:i.blackLevel,dng_whitelevel:i.whiteLevel}}},M=typeof performance<"u"?performance.now():Date.now();return console.info("[Sony cRAW HQ] decode timings",{width:i.width,height:i.height,pyodideReadyMs:Math.round(s-e),lutLoadMs:Math.round(u-s),llvc3DecodeMs:Math.round(c-u),copyMs:Math.round(p-c),totalMs:Math.round(M-e)}),{rawImageData:{data:f,width:i.width,height:i.height,bayerPattern:i.cfaPattern,blackLevels:i.blackLevel,whiteLevel:i.whiteLevel,metadata:b,isThreePlane:!1,isXTrans:!1},info:i}}async function Ws(n,t){return js(n)}var ni=typeof self<"u"?self:global;const Pn=typeof navigator<"u",Hs=Pn&&typeof HTMLImageElement>"u",Qn=!(typeof global>"u"||typeof process>"u"||!process.versions||!process.versions.node),ii=ni.Buffer,On=ni.BigInt,ri=!!ii,Qs=n=>n;function qn(n,t=Qs){if(Qn)try{return typeof require=="function"?Promise.resolve(t(require(n))):import(n).then(t)}catch{console.warn(`Couldn't load ${n}`)}}let zi=ni.fetch;const qs=n=>zi=n;if(!ni.fetch){const n=qn("http",i=>i),t=qn("https",i=>i),e=(i,{headers:r}={})=>new Promise(async(s,o)=>{let{port:a,hostname:l,pathname:u,protocol:d,search:c}=new URL(i);const g={method:"GET",hostname:l,path:encodeURI(u)+c,headers:r};a!==""&&(g.port=Number(a));const h=(d==="https:"?await t:await n).request(g,f=>{if(f.statusCode===301||f.statusCode===302){let p=new URL(f.headers.location,i).toString();return e(p,{headers:r}).then(s).catch(o)}s({status:f.statusCode,arrayBuffer:()=>new Promise(p=>{let m=[];f.on("data",y=>m.push(y)),f.on("end",()=>p(Buffer.concat(m)))})})});h.on("error",o),h.end()});qs(e)}function nt(n,t,e){return t in n?Object.defineProperty(n,t,{value:e,enumerable:!0,configurable:!0,writable:!0}):n[t]=e,n}const Kn=n=>jr(n)?void 0:n,Ks=n=>n!==void 0;function jr(n){return n===void 0||(n instanceof Map?n.size===0:Object.values(n).filter(Ks).length===0)}function kt(n){let t=new Error(n);throw delete t.stack,t}function qe(n){return(n=function(t){for(;t.endsWith("\0");)t=t.slice(0,-1);return t}(n).trim())===""?void 0:n}function Li(n){let t=function(e){let i=0;return e.ifd0.enabled&&(i+=1024),e.exif.enabled&&(i+=2048),e.makerNote&&(i+=2048),e.userComment&&(i+=1024),e.gps.enabled&&(i+=512),e.interop.enabled&&(i+=100),e.ifd1.enabled&&(i+=1024),i+2048}(n);return n.jfif.enabled&&(t+=50),n.xmp.enabled&&(t+=2e4),n.iptc.enabled&&(t+=14e3),n.icc.enabled&&(t+=6e3),t}const Di=n=>String.fromCharCode.apply(null,n),pr=typeof TextDecoder<"u"?new TextDecoder("utf-8"):void 0;function Wr(n){return pr?pr.decode(n):ri?Buffer.from(n).toString("utf8"):decodeURIComponent(escape(Di(n)))}class jt{static from(t,e){return t instanceof this&&t.le===e?t:new jt(t,void 0,void 0,e)}constructor(t,e=0,i,r){if(typeof r=="boolean"&&(this.le=r),Array.isArray(t)&&(t=new Uint8Array(t)),t===0)this.byteOffset=0,this.byteLength=0;else if(t instanceof ArrayBuffer){i===void 0&&(i=t.byteLength-e);let s=new DataView(t,e,i);this._swapDataView(s)}else if(t instanceof Uint8Array||t instanceof DataView||t instanceof jt){i===void 0&&(i=t.byteLength-e),(e+=t.byteOffset)+i>t.byteOffset+t.byteLength&&kt("Creating view outside of available memory in ArrayBuffer");let s=new DataView(t.buffer,e,i);this._swapDataView(s)}else if(typeof t=="number"){let s=new DataView(new ArrayBuffer(t));this._swapDataView(s)}else kt("Invalid input argument for BufferView: "+t)}_swapArrayBuffer(t){this._swapDataView(new DataView(t))}_swapBuffer(t){this._swapDataView(new DataView(t.buffer,t.byteOffset,t.byteLength))}_swapDataView(t){this.dataView=t,this.buffer=t.buffer,this.byteOffset=t.byteOffset,this.byteLength=t.byteLength}_lengthToEnd(t){return this.byteLength-t}set(t,e,i=jt){return t instanceof DataView||t instanceof jt?t=new Uint8Array(t.buffer,t.byteOffset,t.byteLength):t instanceof ArrayBuffer&&(t=new Uint8Array(t)),t instanceof Uint8Array||kt("BufferView.set(): Invalid data argument."),this.toUint8().set(t,e),new i(this,e,t.byteLength)}subarray(t,e){return e=e||this._lengthToEnd(t),new jt(this,t,e)}toUint8(){return new Uint8Array(this.buffer,this.byteOffset,this.byteLength)}getUint8Array(t,e){return new Uint8Array(this.buffer,this.byteOffset+t,e)}getString(t=0,e=this.byteLength){return Wr(this.getUint8Array(t,e))}getLatin1String(t=0,e=this.byteLength){let i=this.getUint8Array(t,e);return Di(i)}getUnicodeString(t=0,e=this.byteLength){const i=[];for(let r=0;r<e&&t+r<this.byteLength;r+=2)i.push(this.getUint16(t+r));return Di(i)}getInt8(t){return this.dataView.getInt8(t)}getUint8(t){return this.dataView.getUint8(t)}getInt16(t,e=this.le){return this.dataView.getInt16(t,e)}getInt32(t,e=this.le){return this.dataView.getInt32(t,e)}getUint16(t,e=this.le){return this.dataView.getUint16(t,e)}getUint32(t,e=this.le){return this.dataView.getUint32(t,e)}getFloat32(t,e=this.le){return this.dataView.getFloat32(t,e)}getFloat64(t,e=this.le){return this.dataView.getFloat64(t,e)}getFloat(t,e=this.le){return this.dataView.getFloat32(t,e)}getDouble(t,e=this.le){return this.dataView.getFloat64(t,e)}getUintBytes(t,e,i){switch(e){case 1:return this.getUint8(t,i);case 2:return this.getUint16(t,i);case 4:return this.getUint32(t,i);case 8:return this.getUint64&&this.getUint64(t,i)}}getUint(t,e,i){switch(e){case 8:return this.getUint8(t,i);case 16:return this.getUint16(t,i);case 32:return this.getUint32(t,i);case 64:return this.getUint64&&this.getUint64(t,i)}}toString(t){return this.dataView.toString(t,this.constructor.name)}ensureChunk(){}}function Bi(n,t){kt(`${n} '${t}' was not loaded, try using full build of exifr.`)}class Vi extends Map{constructor(t){super(),this.kind=t}get(t,e){return this.has(t)||Bi(this.kind,t),e&&(t in e||function(i,r){kt(`Unknown ${i} '${r}'.`)}(this.kind,t),e[t].enabled||Bi(this.kind,t)),super.get(t)}keyList(){return Array.from(this.keys())}}var ye=new Vi("file parser"),Ct=new Vi("segment parser"),be=new Vi("file reader");function $s(n,t){return typeof n=="string"?mr(n,t):Pn&&!Hs&&n instanceof HTMLImageElement?mr(n.src,t):n instanceof Uint8Array||n instanceof ArrayBuffer||n instanceof DataView?new jt(n):Pn&&n instanceof Blob?Ui(n,t,"blob",rn):void kt("Invalid input argument")}function mr(n,t){return(e=n).startsWith("data:")||e.length>1e4?Oi(n,t,"base64"):Qn&&n.includes("://")?Ui(n,t,"url",nn):Qn?Oi(n,t,"fs"):Pn?Ui(n,t,"url",nn):void kt("Invalid input argument");var e}async function Ui(n,t,e,i){return be.has(e)?Oi(n,t,e):i?async function(r,s){let o=await s(r);return new jt(o)}(n,i):void kt(`Parser ${e} is not loaded`)}async function Oi(n,t,e){let i=new(be.get(e))(n,t);return await i.read(),i}const nn=n=>zi(n).then(t=>t.arrayBuffer()),rn=n=>new Promise((t,e)=>{let i=new FileReader;i.onloadend=()=>t(i.result||new ArrayBuffer),i.onerror=e,i.readAsArrayBuffer(n)});class Js extends Map{get tagKeys(){return this.allKeys||(this.allKeys=Array.from(this.keys())),this.allKeys}get tagValues(){return this.allValues||(this.allValues=Array.from(this.values())),this.allValues}}function Pt(n,t,e){let i=new Js;for(let[r,s]of e)i.set(r,s);if(Array.isArray(t))for(let r of t)n.set(r,i);else n.set(t,i);return i}function sn(n,t,e){let i,r=n.get(t);for(i of e)r.set(i[0],i[1])}const Nt=new Map,oe=new Map,Le=new Map,Ae=["chunked","firstChunkSize","firstChunkSizeNode","firstChunkSizeBrowser","chunkSize","chunkLimit"],pn=["jfif","xmp","icc","iptc","ihdr"],an=["tiff",...pn],St=["ifd0","ifd1","exif","gps","interop"],Ie=[...an,...St],Ne=["makerNote","userComment"],mn=["translateKeys","translateValues","reviveValues","multiSegment"],Re=[...mn,"sanitize","mergeOutput","silentErrors"];class Hr{get translate(){return this.translateKeys||this.translateValues||this.reviveValues}}class _n extends Hr{get needed(){return this.enabled||this.deps.size>0}constructor(t,e,i,r){if(super(),nt(this,"enabled",!1),nt(this,"skip",new Set),nt(this,"pick",new Set),nt(this,"deps",new Set),nt(this,"translateKeys",!1),nt(this,"translateValues",!1),nt(this,"reviveValues",!1),this.key=t,this.enabled=e,this.parse=this.enabled,this.applyInheritables(r),this.canBeFiltered=St.includes(t),this.canBeFiltered&&(this.dict=Nt.get(t)),i!==void 0)if(Array.isArray(i))this.parse=this.enabled=!0,this.canBeFiltered&&i.length>0&&this.translateTagSet(i,this.pick);else if(typeof i=="object"){if(this.enabled=!0,this.parse=i.parse!==!1,this.canBeFiltered){let{pick:s,skip:o}=i;s&&s.length>0&&this.translateTagSet(s,this.pick),o&&o.length>0&&this.translateTagSet(o,this.skip)}this.applyInheritables(i)}else i===!0||i===!1?this.parse=this.enabled=i:kt(`Invalid options argument: ${i}`)}applyInheritables(t){let e,i;for(e of mn)i=t[e],i!==void 0&&(this[e]=i)}translateTagSet(t,e){if(this.dict){let i,r,{tagKeys:s,tagValues:o}=this.dict;for(i of t)typeof i=="string"?(r=o.indexOf(i),r===-1&&(r=s.indexOf(Number(i))),r!==-1&&e.add(Number(s[r]))):e.add(i)}else for(let i of t)e.add(i)}finalizeFilters(){!this.enabled&&this.deps.size>0?(this.enabled=!0,$n(this.pick,this.deps)):this.enabled&&this.pick.size>0&&$n(this.pick,this.deps)}}var Ot={jfif:!1,tiff:!0,xmp:!1,icc:!1,iptc:!1,ifd0:!0,ifd1:!1,exif:!0,gps:!0,interop:!1,ihdr:void 0,makerNote:!1,userComment:!1,multiSegment:!1,skip:[],pick:[],translateKeys:!0,translateValues:!0,reviveValues:!0,sanitize:!0,mergeOutput:!0,silentErrors:!0,chunked:!0,firstChunkSize:void 0,firstChunkSizeNode:512,firstChunkSizeBrowser:65536,chunkSize:65536,chunkLimit:5},gr=new Map;class on extends Hr{static useCached(t){let e=gr.get(t);return e!==void 0||(e=new this(t),gr.set(t,e)),e}constructor(t){super(),t===!0?this.setupFromTrue():t===void 0?this.setupFromUndefined():Array.isArray(t)?this.setupFromArray(t):typeof t=="object"?this.setupFromObject(t):kt(`Invalid options argument ${t}`),this.firstChunkSize===void 0&&(this.firstChunkSize=Pn?this.firstChunkSizeBrowser:this.firstChunkSizeNode),this.mergeOutput&&(this.ifd1.enabled=!1),this.filterNestedSegmentTags(),this.traverseTiffDependencyTree(),this.checkLoadedPlugins()}setupFromUndefined(){let t;for(t of Ae)this[t]=Ot[t];for(t of Re)this[t]=Ot[t];for(t of Ne)this[t]=Ot[t];for(t of Ie)this[t]=new _n(t,Ot[t],void 0,this)}setupFromTrue(){let t;for(t of Ae)this[t]=Ot[t];for(t of Re)this[t]=Ot[t];for(t of Ne)this[t]=!0;for(t of Ie)this[t]=new _n(t,!0,void 0,this)}setupFromArray(t){let e;for(e of Ae)this[e]=Ot[e];for(e of Re)this[e]=Ot[e];for(e of Ne)this[e]=Ot[e];for(e of Ie)this[e]=new _n(e,!1,void 0,this);this.setupGlobalFilters(t,void 0,St)}setupFromObject(t){let e;for(e of(St.ifd0=St.ifd0||St.image,St.ifd1=St.ifd1||St.thumbnail,Object.assign(this,t),Ae))this[e]=Ci(t[e],Ot[e]);for(e of Re)this[e]=Ci(t[e],Ot[e]);for(e of Ne)this[e]=Ci(t[e],Ot[e]);for(e of an)this[e]=new _n(e,Ot[e],t[e],this);for(e of St)this[e]=new _n(e,Ot[e],t[e],this.tiff);this.setupGlobalFilters(t.pick,t.skip,St,Ie),t.tiff===!0?this.batchEnableWithBool(St,!0):t.tiff===!1?this.batchEnableWithUserValue(St,t):Array.isArray(t.tiff)?this.setupGlobalFilters(t.tiff,void 0,St):typeof t.tiff=="object"&&this.setupGlobalFilters(t.tiff.pick,t.tiff.skip,St)}batchEnableWithBool(t,e){for(let i of t)this[i].enabled=e}batchEnableWithUserValue(t,e){for(let i of t){let r=e[i];this[i].enabled=r!==!1&&r!==void 0}}setupGlobalFilters(t,e,i,r=i){if(t&&t.length){for(let o of r)this[o].enabled=!1;let s=yr(t,i);for(let[o,a]of s)$n(this[o].pick,a),this[o].enabled=!0}else if(e&&e.length){let s=yr(e,i);for(let[o,a]of s)$n(this[o].skip,a)}}filterNestedSegmentTags(){let{ifd0:t,exif:e,xmp:i,iptc:r,icc:s}=this;this.makerNote?e.deps.add(37500):e.skip.add(37500),this.userComment?e.deps.add(37510):e.skip.add(37510),i.enabled||t.skip.add(700),r.enabled||t.skip.add(33723),s.enabled||t.skip.add(34675)}traverseTiffDependencyTree(){let{ifd0:t,exif:e,gps:i,interop:r}=this;r.needed&&(e.deps.add(40965),t.deps.add(40965)),e.needed&&t.deps.add(34665),i.needed&&t.deps.add(34853),this.tiff.enabled=St.some(s=>this[s].enabled===!0)||this.makerNote||this.userComment;for(let s of St)this[s].finalizeFilters()}get onlyTiff(){return!pn.map(t=>this[t].enabled).some(t=>t===!0)&&this.tiff.enabled}checkLoadedPlugins(){for(let t of an)this[t].enabled&&!Ct.has(t)&&Bi("segment parser",t)}}function yr(n,t){let e,i,r,s,o=[];for(r of t){for(s of(e=Nt.get(r),i=[],e))(n.includes(s[0])||n.includes(s[1]))&&i.push(s[0]);i.length&&o.push([r,i])}return o}function Ci(n,t){return n!==void 0?n:t!==void 0?t:void 0}function $n(n,t){for(let e of t)n.add(e)}nt(on,"default",Ot);class De{constructor(t){nt(this,"parsers",{}),nt(this,"output",{}),nt(this,"errors",[]),nt(this,"pushToErrors",e=>this.errors.push(e)),this.options=on.useCached(t)}async read(t){this.file=await $s(t,this.options)}setup(){if(this.fileParser)return;let{file:t}=this,e=t.getUint16(0);for(let[i,r]of ye)if(r.canHandle(t,e))return this.fileParser=new r(this.options,this.file,this.parsers),t[i]=!0;this.file.close&&this.file.close(),kt("Unknown file format")}async parse(){let{output:t,errors:e}=this;return this.setup(),this.options.silentErrors?(await this.executeParsers().catch(this.pushToErrors),e.push(...this.fileParser.errors)):await this.executeParsers(),this.file.close&&this.file.close(),this.options.silentErrors&&e.length>0&&(t.errors=e),Kn(t)}async executeParsers(){let{output:t}=this;await this.fileParser.parse();let e=Object.values(this.parsers).map(async i=>{let r=await i.parse();i.assignToOutput(t,r)});this.options.silentErrors&&(e=e.map(i=>i.catch(this.pushToErrors))),await Promise.all(e)}async extractThumbnail(){this.setup();let{options:t,file:e}=this,i=Ct.get("tiff",t);var r;if(e.tiff?r={start:0,type:"tiff"}:e.jpeg&&(r=await this.fileParser.getOrFindSegment("tiff")),r===void 0)return;let s=await this.fileParser.ensureSegmentChunk(r),o=this.parsers.tiff=new i(s,t,e),a=await o.extractThumbnail();return e.close&&e.close(),a}}async function si(n,t){let e=new De(t);return await e.read(n),e.parse()}var Zs=Object.freeze({__proto__:null,parse:si,Exifr:De,fileParsers:ye,segmentParsers:Ct,fileReaders:be,tagKeys:Nt,tagValues:oe,tagRevivers:Le,createDictionary:Pt,extendDictionary:sn,fetchUrlAsArrayBuffer:nn,readBlobAsArrayBuffer:rn,chunkedProps:Ae,otherSegments:pn,segments:an,tiffBlocks:St,segmentsAndBlocks:Ie,tiffExtractables:Ne,inheritables:mn,allFormatters:Re,Options:on});class ai{constructor(t,e,i){nt(this,"errors",[]),nt(this,"ensureSegmentChunk",async r=>{let s=r.start,o=r.size||65536;if(this.file.chunked)if(this.file.available(s,o))r.chunk=this.file.subarray(s,o);else try{r.chunk=await this.file.readChunk(s,o)}catch(a){kt(`Couldn't read segment: ${JSON.stringify(r)}. ${a.message}`)}else this.file.byteLength>s+o?r.chunk=this.file.subarray(s,o):r.size===void 0?r.chunk=this.file.subarray(s):kt("Segment unreachable: "+JSON.stringify(r));return r.chunk}),this.extendOptions&&this.extendOptions(t),this.options=t,this.file=e,this.parsers=i}injectSegment(t,e){this.options[t].enabled&&this.createParser(t,e)}createParser(t,e){let i=new(Ct.get(t))(e,this.options,this.file);return this.parsers[t]=i}createParsers(t){for(let e of t){let{type:i,chunk:r}=e,s=this.options[i];if(s&&s.enabled){let o=this.parsers[i];o&&o.append||o||this.createParser(i,r)}}}async readSegments(t){let e=t.map(this.ensureSegmentChunk);await Promise.all(e)}}class ae{static findPosition(t,e){let i=t.getUint16(e+2)+2,r=typeof this.headerLength=="function"?this.headerLength(t,e,i):this.headerLength,s=e+r,o=i-r;return{offset:e,length:i,headerLength:r,start:s,size:o,end:s+o}}static parse(t,e={}){return new this(t,new on({[this.type]:e}),t).parse()}normalizeInput(t){return t instanceof jt?t:new jt(t)}constructor(t,e={},i){nt(this,"errors",[]),nt(this,"raw",new Map),nt(this,"handleError",r=>{if(!this.options.silentErrors)throw r;this.errors.push(r.message)}),this.chunk=this.normalizeInput(t),this.file=i,this.type=this.constructor.type,this.globalOptions=this.options=e,this.localOptions=e[this.type],this.canTranslate=this.localOptions&&this.localOptions.translate}translate(){this.canTranslate&&(this.translated=this.translateBlock(this.raw,this.type))}get output(){return this.translated?this.translated:this.raw?Object.fromEntries(this.raw):void 0}translateBlock(t,e){let i=Le.get(e),r=oe.get(e),s=Nt.get(e),o=this.options[e],a=o.reviveValues&&!!i,l=o.translateValues&&!!r,u=o.translateKeys&&!!s,d={};for(let[c,g]of t)a&&i.has(c)?g=i.get(c)(g):l&&r.has(c)&&(g=this.translateValue(g,r.get(c))),u&&s.has(c)&&(c=s.get(c)||c),d[c]=g;return d}translateValue(t,e){return e[t]||e.DEFAULT||t}assignToOutput(t,e){this.assignObjectToOutput(t,this.constructor.type,e)}assignObjectToOutput(t,e,i){if(this.globalOptions.mergeOutput)return Object.assign(t,i);t[e]?Object.assign(t[e],i):t[e]=i}}nt(ae,"headerLength",4),nt(ae,"type",void 0),nt(ae,"multiSegment",!1),nt(ae,"canHandle",()=>!1);function ta(n){return n===192||n===194||n===196||n===219||n===221||n===218||n===254}function ea(n){return n>=224&&n<=239}function na(n,t,e){for(let[i,r]of Ct)if(r.canHandle(n,t,e))return i}class xr extends ai{constructor(...t){super(...t),nt(this,"appSegments",[]),nt(this,"jpegSegments",[]),nt(this,"unknownSegments",[])}static canHandle(t,e){return e===65496}async parse(){await this.findAppSegments(),await this.readSegments(this.appSegments),this.mergeMultiSegments(),this.createParsers(this.mergedAppSegments||this.appSegments)}setupSegmentFinderArgs(t){t===!0?(this.findAll=!0,this.wanted=new Set(Ct.keyList())):(t=t===void 0?Ct.keyList().filter(e=>this.options[e].enabled):t.filter(e=>this.options[e].enabled&&Ct.has(e)),this.findAll=!1,this.remaining=new Set(t),this.wanted=new Set(t)),this.unfinishedMultiSegment=!1}async findAppSegments(t=0,e){this.setupSegmentFinderArgs(e);let{file:i,findAll:r,wanted:s,remaining:o}=this;if(!r&&this.file.chunked&&(r=Array.from(s).some(a=>{let l=Ct.get(a),u=this.options[a];return l.multiSegment&&u.multiSegment}),r&&await this.file.readWhole()),t=this.findAppSegmentsInRange(t,i.byteLength),!this.options.onlyTiff&&i.chunked){let a=!1;for(;o.size>0&&!a&&(i.canReadNextChunk||this.unfinishedMultiSegment);){let{nextChunkOffset:l}=i,u=this.appSegments.some(d=>!this.file.available(d.offset||d.start,d.length||d.size));if(a=t>l&&!u?!await i.readNextChunk(t):!await i.readNextChunk(l),(t=this.findAppSegmentsInRange(t,i.byteLength))===void 0)return}}}findAppSegmentsInRange(t,e){e-=2;let i,r,s,o,a,l,{file:u,findAll:d,wanted:c,remaining:g,options:h}=this;for(;t<e;t++)if(u.getUint8(t)===255){if(i=u.getUint8(t+1),ea(i)){if(r=u.getUint16(t+2),s=na(u,t,r),s&&c.has(s)&&(o=Ct.get(s),a=o.findPosition(u,t),l=h[s],a.type=s,this.appSegments.push(a),!d&&(o.multiSegment&&l.multiSegment?(this.unfinishedMultiSegment=a.chunkNumber<a.chunkCount,this.unfinishedMultiSegment||g.delete(s)):g.delete(s),g.size===0)))break;h.recordUnknownSegments&&(a=ae.findPosition(u,t),a.marker=i,this.unknownSegments.push(a)),t+=r+1}else if(ta(i)){if(r=u.getUint16(t+2),i===218&&h.stopAfterSos!==!1)return;h.recordJpegSegments&&this.jpegSegments.push({offset:t,length:r,marker:i}),t+=r+1}}return t}mergeMultiSegments(){if(!this.appSegments.some(e=>e.multiSegment))return;let t=function(e,i){let r,s,o,a=new Map;for(let l=0;l<e.length;l++)r=e[l],s=r[i],a.has(s)?o=a.get(s):a.set(s,o=[]),o.push(r);return Array.from(a)}(this.appSegments,"type");this.mergedAppSegments=t.map(([e,i])=>{let r=Ct.get(e,this.options);return r.handleMultiSegments?{type:e,chunk:r.handleMultiSegments(i)}:i[0]})}getSegment(t){return this.appSegments.find(e=>e.type===t)}async getOrFindSegment(t){let e=this.getSegment(t);return e===void 0&&(await this.findAppSegments(0,[t]),e=this.getSegment(t)),e}}nt(xr,"type","jpeg"),ye.set("jpeg",xr);const ia=[void 0,1,1,2,4,8,1,1,2,4,8,4,8,4];class ra extends ae{parseHeader(){var t=this.chunk.getUint16();t===18761?this.le=!0:t===19789&&(this.le=!1),this.chunk.le=this.le,this.headerParsed=!0}parseTags(t,e,i=new Map){let{pick:r,skip:s}=this.options[e];r=new Set(r);let o=r.size>0,a=s.size===0,l=this.chunk.getUint16(t);t+=2;for(let u=0;u<l;u++){let d=this.chunk.getUint16(t);if(o){if(r.has(d)&&(i.set(d,this.parseTag(t,d,e)),r.delete(d),r.size===0))break}else!a&&s.has(d)||i.set(d,this.parseTag(t,d,e));t+=12}return i}parseTag(t,e,i){let{chunk:r}=this,s=r.getUint16(t+2),o=r.getUint32(t+4),a=ia[s];if(a*o<=4?t+=8:t=r.getUint32(t+8),(s<1||s>13)&&kt(`Invalid TIFF value type. block: ${i.toUpperCase()}, tag: ${e.toString(16)}, type: ${s}, offset ${t}`),t>r.byteLength&&kt(`Invalid TIFF value offset. block: ${i.toUpperCase()}, tag: ${e.toString(16)}, type: ${s}, offset ${t} is outside of chunk size ${r.byteLength}`),s===1)return r.getUint8Array(t,o);if(s===2)return qe(r.getString(t,o));if(s===7)return r.getUint8Array(t,o);if(o===1)return this.parseTagValue(s,t);{let l=new(function(d){switch(d){case 1:return Uint8Array;case 3:return Uint16Array;case 4:return Uint32Array;case 5:return Array;case 6:return Int8Array;case 8:return Int16Array;case 9:return Int32Array;case 10:return Array;case 11:return Float32Array;case 12:return Float64Array;default:return Array}}(s))(o),u=a;for(let d=0;d<o;d++)l[d]=this.parseTagValue(s,t),t+=u;return l}}parseTagValue(t,e){let{chunk:i}=this;switch(t){case 1:return i.getUint8(e);case 3:return i.getUint16(e);case 4:return i.getUint32(e);case 5:return i.getUint32(e)/i.getUint32(e+4);case 6:return i.getInt8(e);case 8:return i.getInt16(e);case 9:return i.getInt32(e);case 10:return i.getInt32(e)/i.getInt32(e+4);case 11:return i.getFloat(e);case 12:return i.getDouble(e);case 13:return i.getUint32(e);default:kt(`Invalid tiff type ${t}`)}}}class ki extends ra{static canHandle(t,e){return t.getUint8(e+1)===225&&t.getUint32(e+4)===1165519206&&t.getUint16(e+8)===0}async parse(){this.parseHeader();let{options:t}=this;return t.ifd0.enabled&&await this.parseIfd0Block(),t.exif.enabled&&await this.safeParse("parseExifBlock"),t.gps.enabled&&await this.safeParse("parseGpsBlock"),t.interop.enabled&&await this.safeParse("parseInteropBlock"),t.ifd1.enabled&&await this.safeParse("parseThumbnailBlock"),this.createOutput()}safeParse(t){let e=this[t]();return e.catch!==void 0&&(e=e.catch(this.handleError)),e}findIfd0Offset(){this.ifd0Offset===void 0&&(this.ifd0Offset=this.chunk.getUint32(4))}findIfd1Offset(){if(this.ifd1Offset===void 0){this.findIfd0Offset();let t=this.chunk.getUint16(this.ifd0Offset),e=this.ifd0Offset+2+12*t;this.ifd1Offset=this.chunk.getUint32(e)}}parseBlock(t,e){let i=new Map;return this[e]=i,this.parseTags(t,e,i),i}async parseIfd0Block(){if(this.ifd0)return;let{file:t}=this;this.findIfd0Offset(),this.ifd0Offset<8&&kt("Malformed EXIF data"),!t.chunked&&this.ifd0Offset>t.byteLength&&kt(`IFD0 offset points to outside of file.
this.ifd0Offset: ${this.ifd0Offset}, file.byteLength: ${t.byteLength}`),t.tiff&&await t.ensureChunk(this.ifd0Offset,Li(this.options));let e=this.parseBlock(this.ifd0Offset,"ifd0");return e.size!==0?(this.exifOffset=e.get(34665),this.interopOffset=e.get(40965),this.gpsOffset=e.get(34853),this.xmp=e.get(700),this.iptc=e.get(33723),this.icc=e.get(34675),this.options.sanitize&&(e.delete(34665),e.delete(40965),e.delete(34853),e.delete(700),e.delete(33723),e.delete(34675)),e):void 0}async parseExifBlock(){if(this.exif||(this.ifd0||await this.parseIfd0Block(),this.exifOffset===void 0))return;this.file.tiff&&await this.file.ensureChunk(this.exifOffset,Li(this.options));let t=this.parseBlock(this.exifOffset,"exif");return this.interopOffset||(this.interopOffset=t.get(40965)),this.makerNote=t.get(37500),this.userComment=t.get(37510),this.options.sanitize&&(t.delete(40965),t.delete(37500),t.delete(37510)),this.unpack(t,41728),this.unpack(t,41729),t}unpack(t,e){let i=t.get(e);i&&i.length===1&&t.set(e,i[0])}async parseGpsBlock(){if(this.gps||(this.ifd0||await this.parseIfd0Block(),this.gpsOffset===void 0))return;let t=this.parseBlock(this.gpsOffset,"gps");return t&&t.has(2)&&t.has(4)&&(t.set("latitude",br(...t.get(2),t.get(1))),t.set("longitude",br(...t.get(4),t.get(3)))),t}async parseInteropBlock(){if(!this.interop&&(this.ifd0||await this.parseIfd0Block(),this.interopOffset!==void 0||this.exif||await this.parseExifBlock(),this.interopOffset!==void 0))return this.parseBlock(this.interopOffset,"interop")}async parseThumbnailBlock(t=!1){if(!this.ifd1&&!this.ifd1Parsed&&(!this.options.mergeOutput||t))return this.findIfd1Offset(),this.ifd1Offset>0&&(this.parseBlock(this.ifd1Offset,"ifd1"),this.ifd1Parsed=!0),this.ifd1}async extractThumbnail(){if(this.headerParsed||this.parseHeader(),this.ifd1Parsed||await this.parseThumbnailBlock(!0),this.ifd1===void 0)return;let t=this.ifd1.get(513),e=this.ifd1.get(514);return this.chunk.getUint8Array(t,e)}get image(){return this.ifd0}get thumbnail(){return this.ifd1}createOutput(){let t,e,i,r={};for(e of St)if(t=this[e],!jr(t))if(i=this.canTranslate?this.translateBlock(t,e):Object.fromEntries(t),this.options.mergeOutput){if(e==="ifd1")continue;Object.assign(r,i)}else r[e]=i;return this.makerNote&&(r.makerNote=this.makerNote),this.userComment&&(r.userComment=this.userComment),r}assignToOutput(t,e){if(this.globalOptions.mergeOutput)Object.assign(t,e);else for(let[i,r]of Object.entries(e))this.assignObjectToOutput(t,i,r)}}function br(n,t,e,i){var r=n+t/60+e/3600;return i!=="S"&&i!=="W"||(r*=-1),r}nt(ki,"type","tiff"),nt(ki,"headerLength",10),Ct.set("tiff",ki);var sa=Object.freeze({__proto__:null,default:Zs,Exifr:De,fileParsers:ye,segmentParsers:Ct,fileReaders:be,tagKeys:Nt,tagValues:oe,tagRevivers:Le,createDictionary:Pt,extendDictionary:sn,fetchUrlAsArrayBuffer:nn,readBlobAsArrayBuffer:rn,chunkedProps:Ae,otherSegments:pn,segments:an,tiffBlocks:St,segmentsAndBlocks:Ie,tiffExtractables:Ne,inheritables:mn,allFormatters:Re,Options:on,parse:si});const ji={ifd0:!1,ifd1:!1,exif:!1,gps:!1,interop:!1,sanitize:!1,reviveValues:!0,translateKeys:!1,translateValues:!1,mergeOutput:!1},Wi=Object.assign({},ji,{firstChunkSize:4e4,gps:[1,2,3,4]});async function Qr(n){let t=new De(Wi);await t.read(n);let e=await t.parse();if(e&&e.gps){let{latitude:i,longitude:r}=e.gps;return{latitude:i,longitude:r}}}const Hi=Object.assign({},ji,{tiff:!1,ifd1:!0,mergeOutput:!1});async function qr(n){let t=new De(Hi);await t.read(n);let e=await t.extractThumbnail();return e&&ri?ii.from(e):e}async function Kr(n){let t=await this.thumbnail(n);if(t!==void 0){let e=new Blob([t]);return URL.createObjectURL(e)}}const Qi=Object.assign({},ji,{firstChunkSize:4e4,ifd0:[274]});async function qi(n){let t=new De(Qi);await t.read(n);let e=await t.parse();if(e&&e.ifd0)return e.ifd0[274]}const Ki=Object.freeze({1:{dimensionSwapped:!1,scaleX:1,scaleY:1,deg:0,rad:0},2:{dimensionSwapped:!1,scaleX:-1,scaleY:1,deg:0,rad:0},3:{dimensionSwapped:!1,scaleX:1,scaleY:1,deg:180,rad:180*Math.PI/180},4:{dimensionSwapped:!1,scaleX:-1,scaleY:1,deg:180,rad:180*Math.PI/180},5:{dimensionSwapped:!0,scaleX:1,scaleY:-1,deg:90,rad:90*Math.PI/180},6:{dimensionSwapped:!0,scaleX:1,scaleY:1,deg:90,rad:90*Math.PI/180},7:{dimensionSwapped:!0,scaleX:1,scaleY:-1,deg:270,rad:270*Math.PI/180},8:{dimensionSwapped:!0,scaleX:1,scaleY:1,deg:270,rad:270*Math.PI/180}});let Ke=!0,$e=!0;if(typeof navigator=="object"){let n=navigator.userAgent;if(n.includes("iPad")||n.includes("iPhone")){let t=n.match(/OS (\d+)_(\d+)/);if(t){let[,e,i]=t;Ke=Number(e)+.1*Number(i)<13.4,$e=!1}}else if(n.includes("OS X 10")){let[,t]=n.match(/OS X 10[_.](\d+)/);Ke=$e=Number(t)<15}if(n.includes("Chrome/")){let[,t]=n.match(/Chrome\/(\d+)/);Ke=$e=Number(t)<81}else if(n.includes("Firefox/")){let[,t]=n.match(/Firefox\/(\d+)/);Ke=$e=Number(t)<77}}async function $r(n){let t=await qi(n);return Object.assign({canvas:Ke,css:$e},Ki[t])}class aa extends jt{constructor(...t){super(...t),nt(this,"ranges",new oa),this.byteLength!==0&&this.ranges.add(0,this.byteLength)}_tryExtend(t,e,i){if(t===0&&this.byteLength===0&&i){let r=new DataView(i.buffer||i,i.byteOffset,i.byteLength);this._swapDataView(r)}else{let r=t+e;if(r>this.byteLength){let{dataView:s}=this._extend(r);this._swapDataView(s)}}}_extend(t){let e;e=ri?ii.allocUnsafe(t):new Uint8Array(t);let i=new DataView(e.buffer,e.byteOffset,e.byteLength);return e.set(new Uint8Array(this.buffer,this.byteOffset,this.byteLength),0),{uintView:e,dataView:i}}subarray(t,e,i=!1){return e=e||this._lengthToEnd(t),i&&this._tryExtend(t,e),this.ranges.add(t,e),super.subarray(t,e)}set(t,e,i=!1){i&&this._tryExtend(e,t.byteLength,t);let r=super.set(t,e);return this.ranges.add(e,r.byteLength),r}async ensureChunk(t,e){this.chunked&&(this.ranges.available(t,e)||await this.readChunk(t,e))}available(t,e){return this.ranges.available(t,e)}}class oa{constructor(){nt(this,"list",[])}get length(){return this.list.length}add(t,e,i=0){let r=t+e,s=this.list.filter(o=>Mr(t,o.offset,r)||Mr(t,o.end,r));if(s.length>0){t=Math.min(t,...s.map(a=>a.offset)),r=Math.max(r,...s.map(a=>a.end)),e=r-t;let o=s.shift();o.offset=t,o.length=e,o.end=r,this.list=this.list.filter(a=>!s.includes(a))}else this.list.push({offset:t,length:e,end:r})}available(t,e){let i=t+e;return this.list.some(r=>r.offset<=t&&i<=r.end)}}function Mr(n,t,e){return n<=t&&t<=e}class oi extends aa{constructor(t,e){super(0),nt(this,"chunksRead",0),this.input=t,this.options=e}async readWhole(){this.chunked=!1,await this.readChunk(this.nextChunkOffset)}async readChunked(){this.chunked=!0,await this.readChunk(0,this.options.firstChunkSize)}async readNextChunk(t=this.nextChunkOffset){if(this.fullyRead)return this.chunksRead++,!1;let e=this.options.chunkSize,i=await this.readChunk(t,e);return!!i&&i.byteLength===e}async readChunk(t,e){if(this.chunksRead++,(e=this.safeWrapAddress(t,e))!==0)return this._readChunk(t,e)}safeWrapAddress(t,e){return this.size!==void 0&&t+e>this.size?Math.max(0,this.size-t):e}get nextChunkOffset(){if(this.ranges.list.length!==0)return this.ranges.list[0].length}get canReadNextChunk(){return this.chunksRead<this.options.chunkLimit}get fullyRead(){return this.size!==void 0&&this.nextChunkOffset===this.size}read(){return this.options.chunked?this.readChunked():this.readWhole()}close(){}}be.set("blob",class extends oi{async readWhole(){this.chunked=!1;let n=await rn(this.input);this._swapArrayBuffer(n)}readChunked(){return this.chunked=!0,this.size=this.input.size,super.readChunked()}async _readChunk(n,t){let e=t?n+t:void 0,i=this.input.slice(n,e),r=await rn(i);return this.set(r,n,!0)}});var la=Object.freeze({__proto__:null,default:sa,Exifr:De,fileParsers:ye,segmentParsers:Ct,fileReaders:be,tagKeys:Nt,tagValues:oe,tagRevivers:Le,createDictionary:Pt,extendDictionary:sn,fetchUrlAsArrayBuffer:nn,readBlobAsArrayBuffer:rn,chunkedProps:Ae,otherSegments:pn,segments:an,tiffBlocks:St,segmentsAndBlocks:Ie,tiffExtractables:Ne,inheritables:mn,allFormatters:Re,Options:on,parse:si,gpsOnlyOptions:Wi,gps:Qr,thumbnailOnlyOptions:Hi,thumbnail:qr,thumbnailUrl:Kr,orientationOnlyOptions:Qi,orientation:qi,rotations:Ki,get rotateCanvas(){return Ke},get rotateCss(){return $e},rotation:$r});be.set("url",class extends oi{async readWhole(){this.chunked=!1;let n=await nn(this.input);n instanceof ArrayBuffer?this._swapArrayBuffer(n):n instanceof Uint8Array&&this._swapBuffer(n)}async _readChunk(n,t){let e=t?n+t-1:void 0,i=this.options.httpHeaders||{};(n||e)&&(i.range=`bytes=${[n,e].join("-")}`);let r=await zi(this.input,{headers:i}),s=await r.arrayBuffer(),o=s.byteLength;if(r.status!==416)return o!==t&&(this.size=n+o),this.set(s,n,!0)}});jt.prototype.getUint64=function(n){let t=this.getUint32(n),e=this.getUint32(n+4);return t<1048575?t<<32|e:typeof On!==void 0?(console.warn("Using BigInt because of type 64uint but JS can only handle 53b numbers."),On(t)<<On(32)|On(e)):void kt("Trying to read 64b value but JS can only handle 53b numbers.")};class ca extends ai{parseBoxes(t=0){let e=[];for(;t<this.file.byteLength-4;){let i=this.parseBoxHead(t);if(e.push(i),i.length===0)break;t+=i.length}return e}parseSubBoxes(t){t.boxes=this.parseBoxes(t.start)}findBox(t,e){return t.boxes===void 0&&this.parseSubBoxes(t),t.boxes.find(i=>i.kind===e)}parseBoxHead(t){let e=this.file.getUint32(t),i=this.file.getString(t+4,4),r=t+8;return e===1&&(e=this.file.getUint64(t+8),r+=8),{offset:t,length:e,kind:i,start:r}}parseBoxFullHead(t){if(t.version!==void 0)return;let e=this.file.getUint32(t.start);t.version=e>>24,t.start+=4}}class Jr extends ca{static canHandle(t,e){if(e!==0)return!1;let i=t.getUint16(2);if(i>50)return!1;let r=16,s=[];for(;r<i;)s.push(t.getString(r,4)),r+=4;return s.includes(this.type)}async parse(){let t=this.file.getUint32(0),e=this.parseBoxHead(t);for(;e.kind!=="meta";)t+=e.length,await this.file.ensureChunk(t,16),e=this.parseBoxHead(t);await this.file.ensureChunk(e.offset,e.length),this.parseBoxFullHead(e),this.parseSubBoxes(e),this.options.icc.enabled&&await this.findIcc(e),this.options.tiff.enabled&&await this.findExif(e)}async registerSegment(t,e,i){await this.file.ensureChunk(e,i);let r=this.file.subarray(e,i);this.createParser(t,r)}async findIcc(t){let e=this.findBox(t,"iprp");if(e===void 0)return;let i=this.findBox(e,"ipco");if(i===void 0)return;let r=this.findBox(i,"colr");r!==void 0&&await this.registerSegment("icc",r.offset+12,r.length)}async findExif(t){let e=this.findBox(t,"iinf");if(e===void 0)return;let i=this.findBox(t,"iloc");if(i===void 0)return;let r=this.findExifLocIdInIinf(e),s=this.findExtentInIloc(i,r);if(s===void 0)return;let[o,a]=s;await this.file.ensureChunk(o,a);let l=4+this.file.getUint32(o);o+=l,a-=l,await this.registerSegment("tiff",o,a)}findExifLocIdInIinf(t){this.parseBoxFullHead(t);let e,i,r,s,o=t.start,a=this.file.getUint16(o);for(o+=2;a--;){if(e=this.parseBoxHead(o),this.parseBoxFullHead(e),i=e.start,e.version>=2&&(r=e.version===3?4:2,s=this.file.getString(i+r+2,4),s==="Exif"))return this.file.getUintBytes(i,r);o+=e.length}}get8bits(t){let e=this.file.getUint8(t);return[e>>4,15&e]}findExtentInIloc(t,e){this.parseBoxFullHead(t);let i=t.start,[r,s]=this.get8bits(i++),[o,a]=this.get8bits(i++),l=t.version===2?4:2,u=t.version===1||t.version===2?2:0,d=a+r+s,c=t.version===2?4:2,g=this.file.getUintBytes(i,c);for(i+=c;g--;){let h=this.file.getUintBytes(i,l);i+=l+u+2+o;let f=this.file.getUint16(i);if(i+=2,h===e)return f>1&&console.warn(`ILOC box has more than one extent but we're only processing one
Please create an issue at https://github.com/MikeKovarik/exifr with this file`),[this.file.getUintBytes(i+a,r),this.file.getUintBytes(i+a+r,s)];i+=f*d}}}class Zr extends Jr{}nt(Zr,"type","heic");class _r extends Jr{}nt(_r,"type","avif"),ye.set("heic",Zr),ye.set("avif",_r),Pt(Nt,["ifd0","ifd1"],[[256,"ImageWidth"],[257,"ImageHeight"],[258,"BitsPerSample"],[259,"Compression"],[262,"PhotometricInterpretation"],[270,"ImageDescription"],[271,"Make"],[272,"Model"],[273,"StripOffsets"],[274,"Orientation"],[277,"SamplesPerPixel"],[278,"RowsPerStrip"],[279,"StripByteCounts"],[282,"XResolution"],[283,"YResolution"],[284,"PlanarConfiguration"],[296,"ResolutionUnit"],[301,"TransferFunction"],[305,"Software"],[306,"ModifyDate"],[315,"Artist"],[316,"HostComputer"],[317,"Predictor"],[318,"WhitePoint"],[319,"PrimaryChromaticities"],[513,"ThumbnailOffset"],[514,"ThumbnailLength"],[529,"YCbCrCoefficients"],[530,"YCbCrSubSampling"],[531,"YCbCrPositioning"],[532,"ReferenceBlackWhite"],[700,"ApplicationNotes"],[33432,"Copyright"],[33723,"IPTC"],[34665,"ExifIFD"],[34675,"ICC"],[34853,"GpsIFD"],[330,"SubIFD"],[40965,"InteropIFD"],[40091,"XPTitle"],[40092,"XPComment"],[40093,"XPAuthor"],[40094,"XPKeywords"],[40095,"XPSubject"]]),Pt(Nt,"exif",[[33434,"ExposureTime"],[33437,"FNumber"],[34850,"ExposureProgram"],[34852,"SpectralSensitivity"],[34855,"ISO"],[34858,"TimeZoneOffset"],[34859,"SelfTimerMode"],[34864,"SensitivityType"],[34865,"StandardOutputSensitivity"],[34866,"RecommendedExposureIndex"],[34867,"ISOSpeed"],[34868,"ISOSpeedLatitudeyyy"],[34869,"ISOSpeedLatitudezzz"],[36864,"ExifVersion"],[36867,"DateTimeOriginal"],[36868,"CreateDate"],[36873,"GooglePlusUploadCode"],[36880,"OffsetTime"],[36881,"OffsetTimeOriginal"],[36882,"OffsetTimeDigitized"],[37121,"ComponentsConfiguration"],[37122,"CompressedBitsPerPixel"],[37377,"ShutterSpeedValue"],[37378,"ApertureValue"],[37379,"BrightnessValue"],[37380,"ExposureCompensation"],[37381,"MaxApertureValue"],[37382,"SubjectDistance"],[37383,"MeteringMode"],[37384,"LightSource"],[37385,"Flash"],[37386,"FocalLength"],[37393,"ImageNumber"],[37394,"SecurityClassification"],[37395,"ImageHistory"],[37396,"SubjectArea"],[37500,"MakerNote"],[37510,"UserComment"],[37520,"SubSecTime"],[37521,"SubSecTimeOriginal"],[37522,"SubSecTimeDigitized"],[37888,"AmbientTemperature"],[37889,"Humidity"],[37890,"Pressure"],[37891,"WaterDepth"],[37892,"Acceleration"],[37893,"CameraElevationAngle"],[40960,"FlashpixVersion"],[40961,"ColorSpace"],[40962,"ExifImageWidth"],[40963,"ExifImageHeight"],[40964,"RelatedSoundFile"],[41483,"FlashEnergy"],[41486,"FocalPlaneXResolution"],[41487,"FocalPlaneYResolution"],[41488,"FocalPlaneResolutionUnit"],[41492,"SubjectLocation"],[41493,"ExposureIndex"],[41495,"SensingMethod"],[41728,"FileSource"],[41729,"SceneType"],[41730,"CFAPattern"],[41985,"CustomRendered"],[41986,"ExposureMode"],[41987,"WhiteBalance"],[41988,"DigitalZoomRatio"],[41989,"FocalLengthIn35mmFormat"],[41990,"SceneCaptureType"],[41991,"GainControl"],[41992,"Contrast"],[41993,"Saturation"],[41994,"Sharpness"],[41996,"SubjectDistanceRange"],[42016,"ImageUniqueID"],[42032,"OwnerName"],[42033,"SerialNumber"],[42034,"LensInfo"],[42035,"LensMake"],[42036,"LensModel"],[42037,"LensSerialNumber"],[42080,"CompositeImage"],[42081,"CompositeImageCount"],[42082,"CompositeImageExposureTimes"],[42240,"Gamma"],[59932,"Padding"],[59933,"OffsetSchema"],[65e3,"OwnerName"],[65001,"SerialNumber"],[65002,"Lens"],[65100,"RawFile"],[65101,"Converter"],[65102,"WhiteBalance"],[65105,"Exposure"],[65106,"Shadows"],[65107,"Brightness"],[65108,"Contrast"],[65109,"Saturation"],[65110,"Sharpness"],[65111,"Smoothness"],[65112,"MoireFilter"],[40965,"InteropIFD"]]),Pt(Nt,"gps",[[0,"GPSVersionID"],[1,"GPSLatitudeRef"],[2,"GPSLatitude"],[3,"GPSLongitudeRef"],[4,"GPSLongitude"],[5,"GPSAltitudeRef"],[6,"GPSAltitude"],[7,"GPSTimeStamp"],[8,"GPSSatellites"],[9,"GPSStatus"],[10,"GPSMeasureMode"],[11,"GPSDOP"],[12,"GPSSpeedRef"],[13,"GPSSpeed"],[14,"GPSTrackRef"],[15,"GPSTrack"],[16,"GPSImgDirectionRef"],[17,"GPSImgDirection"],[18,"GPSMapDatum"],[19,"GPSDestLatitudeRef"],[20,"GPSDestLatitude"],[21,"GPSDestLongitudeRef"],[22,"GPSDestLongitude"],[23,"GPSDestBearingRef"],[24,"GPSDestBearing"],[25,"GPSDestDistanceRef"],[26,"GPSDestDistance"],[27,"GPSProcessingMethod"],[28,"GPSAreaInformation"],[29,"GPSDateStamp"],[30,"GPSDifferential"],[31,"GPSHPositioningError"]]),Pt(oe,["ifd0","ifd1"],[[274,{1:"Horizontal (normal)",2:"Mirror horizontal",3:"Rotate 180",4:"Mirror vertical",5:"Mirror horizontal and rotate 270 CW",6:"Rotate 90 CW",7:"Mirror horizontal and rotate 90 CW",8:"Rotate 270 CW"}],[296,{1:"None",2:"inches",3:"cm"}]]);let vn=Pt(oe,"exif",[[34850,{0:"Not defined",1:"Manual",2:"Normal program",3:"Aperture priority",4:"Shutter priority",5:"Creative program",6:"Action program",7:"Portrait mode",8:"Landscape mode"}],[37121,{0:"-",1:"Y",2:"Cb",3:"Cr",4:"R",5:"G",6:"B"}],[37383,{0:"Unknown",1:"Average",2:"CenterWeightedAverage",3:"Spot",4:"MultiSpot",5:"Pattern",6:"Partial",255:"Other"}],[37384,{0:"Unknown",1:"Daylight",2:"Fluorescent",3:"Tungsten (incandescent light)",4:"Flash",9:"Fine weather",10:"Cloudy weather",11:"Shade",12:"Daylight fluorescent (D 5700 - 7100K)",13:"Day white fluorescent (N 4600 - 5400K)",14:"Cool white fluorescent (W 3900 - 4500K)",15:"White fluorescent (WW 3200 - 3700K)",17:"Standard light A",18:"Standard light B",19:"Standard light C",20:"D55",21:"D65",22:"D75",23:"D50",24:"ISO studio tungsten",255:"Other"}],[37385,{0:"Flash did not fire",1:"Flash fired",5:"Strobe return light not detected",7:"Strobe return light detected",9:"Flash fired, compulsory flash mode",13:"Flash fired, compulsory flash mode, return light not detected",15:"Flash fired, compulsory flash mode, return light detected",16:"Flash did not fire, compulsory flash mode",24:"Flash did not fire, auto mode",25:"Flash fired, auto mode",29:"Flash fired, auto mode, return light not detected",31:"Flash fired, auto mode, return light detected",32:"No flash function",65:"Flash fired, red-eye reduction mode",69:"Flash fired, red-eye reduction mode, return light not detected",71:"Flash fired, red-eye reduction mode, return light detected",73:"Flash fired, compulsory flash mode, red-eye reduction mode",77:"Flash fired, compulsory flash mode, red-eye reduction mode, return light not detected",79:"Flash fired, compulsory flash mode, red-eye reduction mode, return light detected",89:"Flash fired, auto mode, red-eye reduction mode",93:"Flash fired, auto mode, return light not detected, red-eye reduction mode",95:"Flash fired, auto mode, return light detected, red-eye reduction mode"}],[41495,{1:"Not defined",2:"One-chip color area sensor",3:"Two-chip color area sensor",4:"Three-chip color area sensor",5:"Color sequential area sensor",7:"Trilinear sensor",8:"Color sequential linear sensor"}],[41728,{1:"Film Scanner",2:"Reflection Print Scanner",3:"Digital Camera"}],[41729,{1:"Directly photographed"}],[41985,{0:"Normal",1:"Custom",2:"HDR (no original saved)",3:"HDR (original saved)",4:"Original (for HDR)",6:"Panorama",7:"Portrait HDR",8:"Portrait"}],[41986,{0:"Auto",1:"Manual",2:"Auto bracket"}],[41987,{0:"Auto",1:"Manual"}],[41990,{0:"Standard",1:"Landscape",2:"Portrait",3:"Night",4:"Other"}],[41991,{0:"None",1:"Low gain up",2:"High gain up",3:"Low gain down",4:"High gain down"}],[41996,{0:"Unknown",1:"Macro",2:"Close",3:"Distant"}],[42080,{0:"Unknown",1:"Not a Composite Image",2:"General Composite Image",3:"Composite Image Captured While Shooting"}]]);const wr={1:"No absolute unit of measurement",2:"Inch",3:"Centimeter"};vn.set(37392,wr),vn.set(41488,wr);const Fi={0:"Normal",1:"Low",2:"High"};function Sr(n){return typeof n=="object"&&n.length!==void 0?n[0]:n}function vr(n){let t=Array.from(n).slice(1);return t[1]>15&&(t=t.map(e=>String.fromCharCode(e))),t[2]!=="0"&&t[2]!==0||t.pop(),t.join(".")}function Ti(n){if(typeof n=="string"){var[t,e,i,r,s,o]=n.trim().split(/[-: ]/g).map(Number),a=new Date(t,e-1,i);return Number.isNaN(r)||Number.isNaN(s)||Number.isNaN(o)||(a.setHours(r),a.setMinutes(s),a.setSeconds(o)),Number.isNaN(+a)?n:a}}function wn(n){if(typeof n=="string")return n;let t=[];if(n[1]===0&&n[n.length-1]===0)for(let e=0;e<n.length;e+=2)t.push(Pr(n[e+1],n[e]));else for(let e=0;e<n.length;e+=2)t.push(Pr(n[e],n[e+1]));return qe(String.fromCodePoint(...t))}function Pr(n,t){return n<<8|t}vn.set(41992,Fi),vn.set(41993,Fi),vn.set(41994,Fi),Pt(Le,["ifd0","ifd1"],[[50827,function(n){return typeof n!="string"?Wr(n):n}],[306,Ti],[40091,wn],[40092,wn],[40093,wn],[40094,wn],[40095,wn]]),Pt(Le,"exif",[[40960,vr],[36864,vr],[36867,Ti],[36868,Ti],[40962,Sr],[40963,Sr]]),Pt(Le,"gps",[[0,n=>Array.from(n).join(".")],[7,n=>Array.from(n).join(":")]]);class Ai extends ae{static canHandle(t,e){return t.getUint8(e+1)===225&&t.getUint32(e+4)===1752462448&&t.getString(e+4,20)==="http://ns.adobe.com/"}static headerLength(t,e){return t.getString(e+4,34)==="http://ns.adobe.com/xmp/extension/"?79:33}static findPosition(t,e){let i=super.findPosition(t,e);return i.multiSegment=i.extended=i.headerLength===79,i.multiSegment?(i.chunkCount=t.getUint8(e+72),i.chunkNumber=t.getUint8(e+76),t.getUint8(e+77)!==0&&i.chunkNumber++):(i.chunkCount=1/0,i.chunkNumber=-1),i}static handleMultiSegments(t){return t.map(e=>e.chunk.getString()).join("")}normalizeInput(t){return typeof t=="string"?t:jt.from(t).getString()}parse(t=this.chunk){if(!this.localOptions.parse)return t;t=function(s){let o={},a={};for(let l of is)o[l]=[],a[l]=0;return s.replace(fa,(l,u,d)=>{if(u==="<"){let c=++a[d];return o[d].push(c),`${l}#${c}`}return`${l}#${o[d].pop()}`})}(t);let e=un.findAll(t,"rdf","Description");e.length===0&&e.push(new un("rdf","Description",void 0,t));let i,r={};for(let s of e)for(let o of s.properties)i=da(o.ns,r),ts(o,i);return function(s){let o;for(let a in s)o=s[a]=Kn(s[a]),o===void 0&&delete s[a];return Kn(s)}(r)}assignToOutput(t,e){if(this.localOptions.parse)for(let[i,r]of Object.entries(e))switch(i){case"tiff":this.assignObjectToOutput(t,"ifd0",r);break;case"exif":this.assignObjectToOutput(t,"exif",r);break;case"xmlns":break;default:this.assignObjectToOutput(t,i,r)}else t.xmp=e}}nt(Ai,"type","xmp"),nt(Ai,"multiSegment",!0),Ct.set("xmp",Ai);class Jn{static findAll(t){return es(t,/([a-zA-Z0-9-]+):([a-zA-Z0-9-]+)=("[^"]*"|'[^']*')/gm).map(Jn.unpackMatch)}static unpackMatch(t){let e=t[1],i=t[2],r=t[3].slice(1,-1);return r=ns(r),new Jn(e,i,r)}constructor(t,e,i){this.ns=t,this.name=e,this.value=i}serialize(){return this.value}}class un{static findAll(t,e,i){if(e!==void 0||i!==void 0){e=e||"[\\w\\d-]+",i=i||"[\\w\\d-]+";var r=new RegExp(`<(${e}):(${i})(#\\d+)?((\\s+?[\\w\\d-:]+=("[^"]*"|'[^']*'))*\\s*)(\\/>|>([\\s\\S]*?)<\\/\\1:\\2\\3>)`,"gm")}else r=/<([\w\d-]+):([\w\d-]+)(#\d+)?((\s+?[\w\d-:]+=("[^"]*"|'[^']*'))*\s*)(\/>|>([\s\S]*?)<\/\1:\2\3>)/gm;return es(t,r).map(un.unpackMatch)}static unpackMatch(t){let e=t[1],i=t[2],r=t[4],s=t[8];return new un(e,i,r,s)}constructor(t,e,i,r){this.ns=t,this.name=e,this.attrString=i,this.innerXml=r,this.attrs=Jn.findAll(i),this.children=un.findAll(r),this.value=this.children.length===0?ns(r):void 0,this.properties=[...this.attrs,...this.children]}get isPrimitive(){return this.value!==void 0&&this.attrs.length===0&&this.children.length===0}get isListContainer(){return this.children.length===1&&this.children[0].isList}get isList(){let{ns:t,name:e}=this;return t==="rdf"&&(e==="Seq"||e==="Bag"||e==="Alt")}get isListItem(){return this.ns==="rdf"&&this.name==="li"}serialize(){if(this.properties.length===0&&this.value===void 0)return;if(this.isPrimitive)return this.value;if(this.isListContainer)return this.children[0].serialize();if(this.isList)return ha(this.children.map(ua));if(this.isListItem&&this.children.length===1&&this.attrs.length===0)return this.children[0].serialize();let t={};for(let e of this.properties)ts(e,t);return this.value!==void 0&&(t.value=this.value),Kn(t)}}function ts(n,t){let e=n.serialize();e!==void 0&&(t[n.name]=e)}var ua=n=>n.serialize(),ha=n=>n.length===1?n[0]:n,da=(n,t)=>t[n]?t[n]:t[n]={};function es(n,t){let e,i=[];if(!n)return i;for(;(e=t.exec(n))!==null;)i.push(e);return i}function ns(n){if(function(i){return i==null||i==="null"||i==="undefined"||i===""||i.trim()===""}(n))return;let t=Number(n);if(!Number.isNaN(t))return t;let e=n.toLowerCase();return e==="true"||e!=="false"&&n.trim()}const is=["rdf:li","rdf:Seq","rdf:Bag","rdf:Alt","rdf:Description"],fa=new RegExp(`(<|\\/)(${is.join("|")})`,"g");var pa=Object.freeze({__proto__:null,default:la,Exifr:De,fileParsers:ye,segmentParsers:Ct,fileReaders:be,tagKeys:Nt,tagValues:oe,tagRevivers:Le,createDictionary:Pt,extendDictionary:sn,fetchUrlAsArrayBuffer:nn,readBlobAsArrayBuffer:rn,chunkedProps:Ae,otherSegments:pn,segments:an,tiffBlocks:St,segmentsAndBlocks:Ie,tiffExtractables:Ne,inheritables:mn,allFormatters:Re,Options:on,parse:si,gpsOnlyOptions:Wi,gps:Qr,thumbnailOnlyOptions:Hi,thumbnail:qr,thumbnailUrl:Kr,orientationOnlyOptions:Qi,orientation:qi,rotations:Ki,get rotateCanvas(){return Ke},get rotateCss(){return $e},rotation:$r});let Cr=qn("fs",n=>n.promises);be.set("fs",class extends oi{async readWhole(){this.chunked=!1,this.fs=await Cr;let n=await this.fs.readFile(this.input);this._swapBuffer(n)}async readChunked(){this.chunked=!0,this.fs=await Cr,await this.open(),await this.readChunk(0,this.options.firstChunkSize)}async open(){this.fh===void 0&&(this.fh=await this.fs.open(this.input,"r"),this.size=(await this.fh.stat(this.input)).size)}async _readChunk(n,t){this.fh===void 0&&await this.open(),n+t>this.size&&(t=this.size-n);var e=this.subarray(n,t,!0);return await this.fh.read(e.dataView,0,t,n),e}async close(){if(this.fh){let n=this.fh;this.fh=void 0,await n.close()}}});be.set("base64",class extends oi{constructor(...n){super(...n),this.input=this.input.replace(/^data:([^;]+);base64,/gim,""),this.size=this.input.length/4*3,this.input.endsWith("==")?this.size-=2:this.input.endsWith("=")&&(this.size-=1)}async _readChunk(n,t){let e,i,r=this.input;n===void 0?(n=0,e=0,i=0):(e=4*Math.floor(n/3),i=n-e/4*3),t===void 0&&(t=this.size);let s=n+t,o=e+4*Math.ceil(s/3);r=r.slice(e,o);let a=Math.min(t,this.size-n);if(ri){let l=ii.from(r,"base64").slice(i,i+a);return this.set(l,n,!0)}{let l=this.subarray(n,a,!0),u=atob(r),d=l.toUint8();for(let c=0;c<a;c++)d[c]=u.charCodeAt(i+c);return l}}});class kr extends ai{static canHandle(t,e){return e===18761||e===19789}extendOptions(t){let{ifd0:e,xmp:i,iptc:r,icc:s}=t;i.enabled&&e.deps.add(700),r.enabled&&e.deps.add(33723),s.enabled&&e.deps.add(34675),e.finalizeFilters()}async parse(){let{tiff:t,xmp:e,iptc:i,icc:r}=this.options;if(t.enabled||e.enabled||i.enabled||r.enabled){let s=Math.max(Li(this.options),this.options.chunkSize);await this.file.ensureChunk(0,s),this.createParser("tiff",this.file),this.parsers.tiff.parseHeader(),await this.parsers.tiff.parseIfd0Block(),this.adaptTiffPropAsSegment("xmp"),this.adaptTiffPropAsSegment("iptc"),this.adaptTiffPropAsSegment("icc")}}adaptTiffPropAsSegment(t){if(this.parsers.tiff[t]){let e=this.parsers.tiff[t];this.injectSegment(t,e)}}}nt(kr,"type","tiff"),ye.set("tiff",kr);let ma=qn("zlib");const ga=["ihdr","iccp","text","itxt","exif"];class Fr extends ai{constructor(...t){super(...t),nt(this,"catchError",e=>this.errors.push(e)),nt(this,"metaChunks",[]),nt(this,"unknownChunks",[])}static canHandle(t,e){return e===35152&&t.getUint32(0)===2303741511&&t.getUint32(4)===218765834}async parse(){let{file:t}=this;await this.findPngChunksInRange(8,t.byteLength),await this.readSegments(this.metaChunks),this.findIhdr(),this.parseTextChunks(),await this.findExif().catch(this.catchError),await this.findXmp().catch(this.catchError),await this.findIcc().catch(this.catchError)}async findPngChunksInRange(t,e){let{file:i}=this;for(;t<e;){let r=i.getUint32(t),s=i.getUint32(t+4),o=i.getString(t+4,4).toLowerCase(),a=r+4+4+4,l={type:o,offset:t,length:a,start:t+4+4,size:r,marker:s};ga.includes(o)?this.metaChunks.push(l):this.unknownChunks.push(l),t+=a}}parseTextChunks(){let t=this.metaChunks.filter(e=>e.type==="text");for(let e of t){let[i,r]=this.file.getString(e.start,e.size).split("\0");this.injectKeyValToIhdr(i,r)}}injectKeyValToIhdr(t,e){let i=this.parsers.ihdr;i&&i.raw.set(t,e)}findIhdr(){let t=this.metaChunks.find(e=>e.type==="ihdr");t&&this.options.ihdr.enabled!==!1&&this.createParser("ihdr",t.chunk)}async findExif(){let t=this.metaChunks.find(e=>e.type==="exif");t&&this.injectSegment("tiff",t.chunk)}async findXmp(){let t=this.metaChunks.filter(e=>e.type==="itxt");for(let e of t)e.chunk.getString(0,17)==="XML:com.adobe.xmp"&&this.injectSegment("xmp",e.chunk)}async findIcc(){let t=this.metaChunks.find(a=>a.type==="iccp");if(!t)return;let{chunk:e}=t,i=e.getUint8Array(0,81),r=0;for(;r<80&&i[r]!==0;)r++;let s=r+2,o=e.getString(0,r);if(this.injectKeyValToIhdr("ProfileName",o),Qn){let a=await ma,l=e.getUint8Array(s);l=a.inflateSync(l),this.injectSegment("icc",l)}}}nt(Fr,"type","png"),ye.set("png",Fr),Pt(Nt,"interop",[[1,"InteropIndex"],[2,"InteropVersion"],[4096,"RelatedImageFileFormat"],[4097,"RelatedImageWidth"],[4098,"RelatedImageHeight"]]),sn(Nt,"ifd0",[[11,"ProcessingSoftware"],[254,"SubfileType"],[255,"OldSubfileType"],[263,"Thresholding"],[264,"CellWidth"],[265,"CellLength"],[266,"FillOrder"],[269,"DocumentName"],[280,"MinSampleValue"],[281,"MaxSampleValue"],[285,"PageName"],[286,"XPosition"],[287,"YPosition"],[290,"GrayResponseUnit"],[297,"PageNumber"],[321,"HalftoneHints"],[322,"TileWidth"],[323,"TileLength"],[332,"InkSet"],[337,"TargetPrinter"],[18246,"Rating"],[18249,"RatingPercent"],[33550,"PixelScale"],[34264,"ModelTransform"],[34377,"PhotoshopSettings"],[50706,"DNGVersion"],[50707,"DNGBackwardVersion"],[50708,"UniqueCameraModel"],[50709,"LocalizedCameraModel"],[50736,"DNGLensInfo"],[50739,"ShadowScale"],[50740,"DNGPrivateData"],[33920,"IntergraphMatrix"],[33922,"ModelTiePoint"],[34118,"SEMInfo"],[34735,"GeoTiffDirectory"],[34736,"GeoTiffDoubleParams"],[34737,"GeoTiffAsciiParams"],[50341,"PrintIM"],[50721,"ColorMatrix1"],[50722,"ColorMatrix2"],[50723,"CameraCalibration1"],[50724,"CameraCalibration2"],[50725,"ReductionMatrix1"],[50726,"ReductionMatrix2"],[50727,"AnalogBalance"],[50728,"AsShotNeutral"],[50729,"AsShotWhiteXY"],[50730,"BaselineExposure"],[50731,"BaselineNoise"],[50732,"BaselineSharpness"],[50734,"LinearResponseLimit"],[50735,"CameraSerialNumber"],[50741,"MakerNoteSafety"],[50778,"CalibrationIlluminant1"],[50779,"CalibrationIlluminant2"],[50781,"RawDataUniqueID"],[50827,"OriginalRawFileName"],[50828,"OriginalRawFileData"],[50831,"AsShotICCProfile"],[50832,"AsShotPreProfileMatrix"],[50833,"CurrentICCProfile"],[50834,"CurrentPreProfileMatrix"],[50879,"ColorimetricReference"],[50885,"SRawType"],[50898,"PanasonicTitle"],[50899,"PanasonicTitle2"],[50931,"CameraCalibrationSig"],[50932,"ProfileCalibrationSig"],[50933,"ProfileIFD"],[50934,"AsShotProfileName"],[50936,"ProfileName"],[50937,"ProfileHueSatMapDims"],[50938,"ProfileHueSatMapData1"],[50939,"ProfileHueSatMapData2"],[50940,"ProfileToneCurve"],[50941,"ProfileEmbedPolicy"],[50942,"ProfileCopyright"],[50964,"ForwardMatrix1"],[50965,"ForwardMatrix2"],[50966,"PreviewApplicationName"],[50967,"PreviewApplicationVersion"],[50968,"PreviewSettingsName"],[50969,"PreviewSettingsDigest"],[50970,"PreviewColorSpace"],[50971,"PreviewDateTime"],[50972,"RawImageDigest"],[50973,"OriginalRawFileDigest"],[50981,"ProfileLookTableDims"],[50982,"ProfileLookTableData"],[51043,"TimeCodes"],[51044,"FrameRate"],[51058,"TStop"],[51081,"ReelName"],[51089,"OriginalDefaultFinalSize"],[51090,"OriginalBestQualitySize"],[51091,"OriginalDefaultCropSize"],[51105,"CameraLabel"],[51107,"ProfileHueSatMapEncoding"],[51108,"ProfileLookTableEncoding"],[51109,"BaselineExposureOffset"],[51110,"DefaultBlackRender"],[51111,"NewRawImageDigest"],[51112,"RawToPreviewGain"]]);let Tr=[[273,"StripOffsets"],[279,"StripByteCounts"],[288,"FreeOffsets"],[289,"FreeByteCounts"],[291,"GrayResponseCurve"],[292,"T4Options"],[293,"T6Options"],[300,"ColorResponseUnit"],[320,"ColorMap"],[324,"TileOffsets"],[325,"TileByteCounts"],[326,"BadFaxLines"],[327,"CleanFaxData"],[328,"ConsecutiveBadFaxLines"],[330,"SubIFD"],[333,"InkNames"],[334,"NumberofInks"],[336,"DotRange"],[338,"ExtraSamples"],[339,"SampleFormat"],[340,"SMinSampleValue"],[341,"SMaxSampleValue"],[342,"TransferRange"],[343,"ClipPath"],[344,"XClipPathUnits"],[345,"YClipPathUnits"],[346,"Indexed"],[347,"JPEGTables"],[351,"OPIProxy"],[400,"GlobalParametersIFD"],[401,"ProfileType"],[402,"FaxProfile"],[403,"CodingMethods"],[404,"VersionYear"],[405,"ModeNumber"],[433,"Decode"],[434,"DefaultImageColor"],[435,"T82Options"],[437,"JPEGTables"],[512,"JPEGProc"],[515,"JPEGRestartInterval"],[517,"JPEGLosslessPredictors"],[518,"JPEGPointTransforms"],[519,"JPEGQTables"],[520,"JPEGDCTables"],[521,"JPEGACTables"],[559,"StripRowCounts"],[999,"USPTOMiscellaneous"],[18247,"XP_DIP_XML"],[18248,"StitchInfo"],[28672,"SonyRawFileType"],[28688,"SonyToneCurve"],[28721,"VignettingCorrection"],[28722,"VignettingCorrParams"],[28724,"ChromaticAberrationCorrection"],[28725,"ChromaticAberrationCorrParams"],[28726,"DistortionCorrection"],[28727,"DistortionCorrParams"],[29895,"SonyCropTopLeft"],[29896,"SonyCropSize"],[32781,"ImageID"],[32931,"WangTag1"],[32932,"WangAnnotation"],[32933,"WangTag3"],[32934,"WangTag4"],[32953,"ImageReferencePoints"],[32954,"RegionXformTackPoint"],[32955,"WarpQuadrilateral"],[32956,"AffineTransformMat"],[32995,"Matteing"],[32996,"DataType"],[32997,"ImageDepth"],[32998,"TileDepth"],[33300,"ImageFullWidth"],[33301,"ImageFullHeight"],[33302,"TextureFormat"],[33303,"WrapModes"],[33304,"FovCot"],[33305,"MatrixWorldToScreen"],[33306,"MatrixWorldToCamera"],[33405,"Model2"],[33421,"CFARepeatPatternDim"],[33422,"CFAPattern2"],[33423,"BatteryLevel"],[33424,"KodakIFD"],[33445,"MDFileTag"],[33446,"MDScalePixel"],[33447,"MDColorTable"],[33448,"MDLabName"],[33449,"MDSampleInfo"],[33450,"MDPrepDate"],[33451,"MDPrepTime"],[33452,"MDFileUnits"],[33589,"AdventScale"],[33590,"AdventRevision"],[33628,"UIC1Tag"],[33629,"UIC2Tag"],[33630,"UIC3Tag"],[33631,"UIC4Tag"],[33918,"IntergraphPacketData"],[33919,"IntergraphFlagRegisters"],[33921,"INGRReserved"],[34016,"Site"],[34017,"ColorSequence"],[34018,"IT8Header"],[34019,"RasterPadding"],[34020,"BitsPerRunLength"],[34021,"BitsPerExtendedRunLength"],[34022,"ColorTable"],[34023,"ImageColorIndicator"],[34024,"BackgroundColorIndicator"],[34025,"ImageColorValue"],[34026,"BackgroundColorValue"],[34027,"PixelIntensityRange"],[34028,"TransparencyIndicator"],[34029,"ColorCharacterization"],[34030,"HCUsage"],[34031,"TrapIndicator"],[34032,"CMYKEquivalent"],[34152,"AFCP_IPTC"],[34232,"PixelMagicJBIGOptions"],[34263,"JPLCartoIFD"],[34306,"WB_GRGBLevels"],[34310,"LeafData"],[34687,"TIFF_FXExtensions"],[34688,"MultiProfiles"],[34689,"SharedData"],[34690,"T88Options"],[34732,"ImageLayer"],[34750,"JBIGOptions"],[34856,"Opto-ElectricConvFactor"],[34857,"Interlace"],[34908,"FaxRecvParams"],[34909,"FaxSubAddress"],[34910,"FaxRecvTime"],[34929,"FedexEDR"],[34954,"LeafSubIFD"],[37387,"FlashEnergy"],[37388,"SpatialFrequencyResponse"],[37389,"Noise"],[37390,"FocalPlaneXResolution"],[37391,"FocalPlaneYResolution"],[37392,"FocalPlaneResolutionUnit"],[37397,"ExposureIndex"],[37398,"TIFF-EPStandardID"],[37399,"SensingMethod"],[37434,"CIP3DataFile"],[37435,"CIP3Sheet"],[37436,"CIP3Side"],[37439,"StoNits"],[37679,"MSDocumentText"],[37680,"MSPropertySetStorage"],[37681,"MSDocumentTextPosition"],[37724,"ImageSourceData"],[40965,"InteropIFD"],[40976,"SamsungRawPointersOffset"],[40977,"SamsungRawPointersLength"],[41217,"SamsungRawByteOrder"],[41218,"SamsungRawUnknown"],[41484,"SpatialFrequencyResponse"],[41485,"Noise"],[41489,"ImageNumber"],[41490,"SecurityClassification"],[41491,"ImageHistory"],[41494,"TIFF-EPStandardID"],[41995,"DeviceSettingDescription"],[42112,"GDALMetadata"],[42113,"GDALNoData"],[44992,"ExpandSoftware"],[44993,"ExpandLens"],[44994,"ExpandFilm"],[44995,"ExpandFilterLens"],[44996,"ExpandScanner"],[44997,"ExpandFlashLamp"],[46275,"HasselbladRawImage"],[48129,"PixelFormat"],[48130,"Transformation"],[48131,"Uncompressed"],[48132,"ImageType"],[48256,"ImageWidth"],[48257,"ImageHeight"],[48258,"WidthResolution"],[48259,"HeightResolution"],[48320,"ImageOffset"],[48321,"ImageByteCount"],[48322,"AlphaOffset"],[48323,"AlphaByteCount"],[48324,"ImageDataDiscard"],[48325,"AlphaDataDiscard"],[50215,"OceScanjobDesc"],[50216,"OceApplicationSelector"],[50217,"OceIDNumber"],[50218,"OceImageLogic"],[50255,"Annotations"],[50459,"HasselbladExif"],[50547,"OriginalFileName"],[50560,"USPTOOriginalContentType"],[50656,"CR2CFAPattern"],[50710,"CFAPlaneColor"],[50711,"CFALayout"],[50712,"LinearizationTable"],[50713,"BlackLevelRepeatDim"],[50714,"BlackLevel"],[50715,"BlackLevelDeltaH"],[50716,"BlackLevelDeltaV"],[50717,"WhiteLevel"],[50718,"DefaultScale"],[50719,"DefaultCropOrigin"],[50720,"DefaultCropSize"],[50733,"BayerGreenSplit"],[50737,"ChromaBlurRadius"],[50738,"AntiAliasStrength"],[50752,"RawImageSegmentation"],[50780,"BestQualityScale"],[50784,"AliasLayerMetadata"],[50829,"ActiveArea"],[50830,"MaskedAreas"],[50935,"NoiseReductionApplied"],[50974,"SubTileBlockSize"],[50975,"RowInterleaveFactor"],[51008,"OpcodeList1"],[51009,"OpcodeList2"],[51022,"OpcodeList3"],[51041,"NoiseProfile"],[51114,"CacheVersion"],[51125,"DefaultUserCrop"],[51157,"NikonNEFInfo"],[65024,"KdcIFD"]];sn(Nt,"ifd0",Tr),sn(Nt,"exif",Tr),Pt(oe,"gps",[[23,{M:"Magnetic North",T:"True North"}],[25,{K:"Kilometers",M:"Miles",N:"Nautical Miles"}]]);class Ii extends ae{static canHandle(t,e){return t.getUint8(e+1)===224&&t.getUint32(e+4)===1246120262&&t.getUint8(e+8)===0}parse(){return this.parseTags(),this.translate(),this.output}parseTags(){this.raw=new Map([[0,this.chunk.getUint16(0)],[2,this.chunk.getUint8(2)],[3,this.chunk.getUint16(3)],[5,this.chunk.getUint16(5)],[7,this.chunk.getUint8(7)],[8,this.chunk.getUint8(8)]])}}nt(Ii,"type","jfif"),nt(Ii,"headerLength",9),Ct.set("jfif",Ii),Pt(Nt,"jfif",[[0,"JFIFVersion"],[2,"ResolutionUnit"],[3,"XResolution"],[5,"YResolution"],[7,"ThumbnailWidth"],[8,"ThumbnailHeight"]]);class Ar extends ae{parse(){return this.parseTags(),this.translate(),this.output}parseTags(){this.raw=new Map([[0,this.chunk.getUint32(0)],[4,this.chunk.getUint32(4)],[8,this.chunk.getUint8(8)],[9,this.chunk.getUint8(9)],[10,this.chunk.getUint8(10)],[11,this.chunk.getUint8(11)],[12,this.chunk.getUint8(12)],...Array.from(this.raw)])}}nt(Ar,"type","ihdr"),Ct.set("ihdr",Ar),Pt(Nt,"ihdr",[[0,"ImageWidth"],[4,"ImageHeight"],[8,"BitDepth"],[9,"ColorType"],[10,"Compression"],[11,"Filter"],[12,"Interlace"]]),Pt(oe,"ihdr",[[9,{0:"Grayscale",2:"RGB",3:"Palette",4:"Grayscale with Alpha",6:"RGB with Alpha",DEFAULT:"Unknown"}],[10,{0:"Deflate/Inflate",DEFAULT:"Unknown"}],[11,{0:"Adaptive",DEFAULT:"Unknown"}],[12,{0:"Noninterlaced",1:"Adam7 Interlace",DEFAULT:"Unknown"}]]);class Hn extends ae{static canHandle(t,e){return t.getUint8(e+1)===226&&t.getUint32(e+4)===1229144927}static findPosition(t,e){let i=super.findPosition(t,e);return i.chunkNumber=t.getUint8(e+16),i.chunkCount=t.getUint8(e+17),i.multiSegment=i.chunkCount>1,i}static handleMultiSegments(t){return function(e){let i=function(r){let s=r[0].constructor,o=0;for(let u of r)o+=u.length;let a=new s(o),l=0;for(let u of r)a.set(u,l),l+=u.length;return a}(e.map(r=>r.chunk.toUint8()));return new jt(i)}(t)}parse(){return this.raw=new Map,this.parseHeader(),this.parseTags(),this.translate(),this.output}parseHeader(){let{raw:t}=this;this.chunk.byteLength<84&&kt("ICC header is too short");for(let[e,i]of Object.entries(ya)){e=parseInt(e,10);let r=i(this.chunk,e);r!=="\0\0\0\0"&&t.set(e,r)}}parseTags(){let t,e,i,r,s,{raw:o}=this,a=this.chunk.getUint32(128),l=132,u=this.chunk.byteLength;for(;a--;){if(t=this.chunk.getString(l,4),e=this.chunk.getUint32(l+4),i=this.chunk.getUint32(l+8),r=this.chunk.getString(e,4),e+i>u)return void console.warn("reached the end of the first ICC chunk. Enable options.tiff.multiSegment to read all ICC segments.");s=this.parseTag(r,e,i),s!==void 0&&s!=="\0\0\0\0"&&o.set(t,s),l+=12}}parseTag(t,e,i){switch(t){case"desc":return this.parseDesc(e);case"mluc":return this.parseMluc(e);case"text":return this.parseText(e,i);case"sig ":return this.parseSig(e)}if(!(e+i>this.chunk.byteLength))return this.chunk.getUint8Array(e,i)}parseDesc(t){let e=this.chunk.getUint32(t+8)-1;return qe(this.chunk.getString(t+12,e))}parseText(t,e){return qe(this.chunk.getString(t+8,e-8))}parseSig(t){return qe(this.chunk.getString(t+8,4))}parseMluc(t){let{chunk:e}=this,i=e.getUint32(t+8),r=e.getUint32(t+12),s=t+16,o=[];for(let a=0;a<i;a++){let l=e.getString(s+0,2),u=e.getString(s+2,2),d=e.getUint32(s+4),c=e.getUint32(s+8)+t,g=qe(e.getUnicodeString(c,d));o.push({lang:l,country:u,text:g}),s+=r}return i===1?o[0].text:o}translateValue(t,e){return typeof t=="string"?e[t]||e[t.toLowerCase()]||t:e[t]||t}}nt(Hn,"type","icc"),nt(Hn,"multiSegment",!0),nt(Hn,"headerLength",18);const ya={4:we,8:function(n,t){return[n.getUint8(t),n.getUint8(t+1)>>4,n.getUint8(t+1)%16].map(e=>e.toString(10)).join(".")},12:we,16:we,20:we,24:function(n,t){const e=n.getUint16(t),i=n.getUint16(t+2)-1,r=n.getUint16(t+4),s=n.getUint16(t+6),o=n.getUint16(t+8),a=n.getUint16(t+10);return new Date(Date.UTC(e,i,r,s,o,a))},36:we,40:we,48:we,52:we,64:(n,t)=>n.getUint32(t),80:we};function we(n,t){return qe(n.getString(t,4))}Ct.set("icc",Hn),Pt(Nt,"icc",[[4,"ProfileCMMType"],[8,"ProfileVersion"],[12,"ProfileClass"],[16,"ColorSpaceData"],[20,"ProfileConnectionSpace"],[24,"ProfileDateTime"],[36,"ProfileFileSignature"],[40,"PrimaryPlatform"],[44,"CMMFlags"],[48,"DeviceManufacturer"],[52,"DeviceModel"],[56,"DeviceAttributes"],[64,"RenderingIntent"],[68,"ConnectionSpaceIlluminant"],[80,"ProfileCreator"],[84,"ProfileID"],["Header","ProfileHeader"],["MS00","WCSProfiles"],["bTRC","BlueTRC"],["bXYZ","BlueMatrixColumn"],["bfd","UCRBG"],["bkpt","MediaBlackPoint"],["calt","CalibrationDateTime"],["chad","ChromaticAdaptation"],["chrm","Chromaticity"],["ciis","ColorimetricIntentImageState"],["clot","ColorantTableOut"],["clro","ColorantOrder"],["clrt","ColorantTable"],["cprt","ProfileCopyright"],["crdi","CRDInfo"],["desc","ProfileDescription"],["devs","DeviceSettings"],["dmdd","DeviceModelDesc"],["dmnd","DeviceMfgDesc"],["dscm","ProfileDescriptionML"],["fpce","FocalPlaneColorimetryEstimates"],["gTRC","GreenTRC"],["gXYZ","GreenMatrixColumn"],["gamt","Gamut"],["kTRC","GrayTRC"],["lumi","Luminance"],["meas","Measurement"],["meta","Metadata"],["mmod","MakeAndModel"],["ncl2","NamedColor2"],["ncol","NamedColor"],["ndin","NativeDisplayInfo"],["pre0","Preview0"],["pre1","Preview1"],["pre2","Preview2"],["ps2i","PS2RenderingIntent"],["ps2s","PostScript2CSA"],["psd0","PostScript2CRD0"],["psd1","PostScript2CRD1"],["psd2","PostScript2CRD2"],["psd3","PostScript2CRD3"],["pseq","ProfileSequenceDesc"],["psid","ProfileSequenceIdentifier"],["psvm","PS2CRDVMSize"],["rTRC","RedTRC"],["rXYZ","RedMatrixColumn"],["resp","OutputResponse"],["rhoc","ReflectionHardcopyOrigColorimetry"],["rig0","PerceptualRenderingIntentGamut"],["rig2","SaturationRenderingIntentGamut"],["rpoc","ReflectionPrintOutputColorimetry"],["sape","SceneAppearanceEstimates"],["scoe","SceneColorimetryEstimates"],["scrd","ScreeningDesc"],["scrn","Screening"],["targ","CharTarget"],["tech","Technology"],["vcgt","VideoCardGamma"],["view","ViewingConditions"],["vued","ViewingCondDesc"],["wtpt","MediaWhitePoint"]]);const Xn={"4d2p":"Erdt Systems",AAMA:"Aamazing Technologies",ACER:"Acer",ACLT:"Acolyte Color Research",ACTI:"Actix Sytems",ADAR:"Adara Technology",ADBE:"Adobe",ADI:"ADI Systems",AGFA:"Agfa Graphics",ALMD:"Alps Electric",ALPS:"Alps Electric",ALWN:"Alwan Color Expertise",AMTI:"Amiable Technologies",AOC:"AOC International",APAG:"Apago",APPL:"Apple Computer",AST:"AST","AT&T":"AT&T",BAEL:"BARBIERI electronic",BRCO:"Barco NV",BRKP:"Breakpoint",BROT:"Brother",BULL:"Bull",BUS:"Bus Computer Systems","C-IT":"C-Itoh",CAMR:"Intel",CANO:"Canon",CARR:"Carroll Touch",CASI:"Casio",CBUS:"Colorbus PL",CEL:"Crossfield",CELx:"Crossfield",CGS:"CGS Publishing Technologies International",CHM:"Rochester Robotics",CIGL:"Colour Imaging Group, London",CITI:"Citizen",CL00:"Candela",CLIQ:"Color IQ",CMCO:"Chromaco",CMiX:"CHROMiX",COLO:"Colorgraphic Communications",COMP:"Compaq",COMp:"Compeq/Focus Technology",CONR:"Conrac Display Products",CORD:"Cordata Technologies",CPQ:"Compaq",CPRO:"ColorPro",CRN:"Cornerstone",CTX:"CTX International",CVIS:"ColorVision",CWC:"Fujitsu Laboratories",DARI:"Darius Technology",DATA:"Dataproducts",DCP:"Dry Creek Photo",DCRC:"Digital Contents Resource Center, Chung-Ang University",DELL:"Dell Computer",DIC:"Dainippon Ink and Chemicals",DICO:"Diconix",DIGI:"Digital","DL&C":"Digital Light & Color",DPLG:"Doppelganger",DS:"Dainippon Screen",DSOL:"DOOSOL",DUPN:"DuPont",EPSO:"Epson",ESKO:"Esko-Graphics",ETRI:"Electronics and Telecommunications Research Institute",EVER:"Everex Systems",EXAC:"ExactCODE",Eizo:"Eizo",FALC:"Falco Data Products",FF:"Fuji Photo Film",FFEI:"FujiFilm Electronic Imaging",FNRD:"Fnord Software",FORA:"Fora",FORE:"Forefront Technology",FP:"Fujitsu",FPA:"WayTech Development",FUJI:"Fujitsu",FX:"Fuji Xerox",GCC:"GCC Technologies",GGSL:"Global Graphics Software",GMB:"Gretagmacbeth",GMG:"GMG",GOLD:"GoldStar Technology",GOOG:"Google",GPRT:"Giantprint",GTMB:"Gretagmacbeth",GVC:"WayTech Development",GW2K:"Sony",HCI:"HCI",HDM:"Heidelberger Druckmaschinen",HERM:"Hermes",HITA:"Hitachi America",HP:"Hewlett-Packard",HTC:"Hitachi",HiTi:"HiTi Digital",IBM:"IBM",IDNT:"Scitex",IEC:"Hewlett-Packard",IIYA:"Iiyama North America",IKEG:"Ikegami Electronics",IMAG:"Image Systems",IMI:"Ingram Micro",INTC:"Intel",INTL:"N/A (INTL)",INTR:"Intra Electronics",IOCO:"Iocomm International Technology",IPS:"InfoPrint Solutions Company",IRIS:"Scitex",ISL:"Ichikawa Soft Laboratory",ITNL:"N/A (ITNL)",IVM:"IVM",IWAT:"Iwatsu Electric",Idnt:"Scitex",Inca:"Inca Digital Printers",Iris:"Scitex",JPEG:"Joint Photographic Experts Group",JSFT:"Jetsoft Development",JVC:"JVC Information Products",KART:"Scitex",KFC:"KFC Computek Components",KLH:"KLH Computers",KMHD:"Konica Minolta",KNCA:"Konica",KODA:"Kodak",KYOC:"Kyocera",Kart:"Scitex",LCAG:"Leica",LCCD:"Leeds Colour",LDAK:"Left Dakota",LEAD:"Leading Technology",LEXM:"Lexmark International",LINK:"Link Computer",LINO:"Linotronic",LITE:"Lite-On",Leaf:"Leaf",Lino:"Linotronic",MAGC:"Mag Computronic",MAGI:"MAG Innovision",MANN:"Mannesmann",MICN:"Micron Technology",MICR:"Microtek",MICV:"Microvitec",MINO:"Minolta",MITS:"Mitsubishi Electronics America",MITs:"Mitsuba",MNLT:"Minolta",MODG:"Modgraph",MONI:"Monitronix",MONS:"Monaco Systems",MORS:"Morse Technology",MOTI:"Motive Systems",MSFT:"Microsoft",MUTO:"MUTOH INDUSTRIES",Mits:"Mitsubishi Electric",NANA:"NANAO",NEC:"NEC",NEXP:"NexPress Solutions",NISS:"Nissei Sangyo America",NKON:"Nikon",NONE:"none",OCE:"Oce Technologies",OCEC:"OceColor",OKI:"Oki",OKID:"Okidata",OKIP:"Okidata",OLIV:"Olivetti",OLYM:"Olympus",ONYX:"Onyx Graphics",OPTI:"Optiquest",PACK:"Packard Bell",PANA:"Matsushita Electric Industrial",PANT:"Pantone",PBN:"Packard Bell",PFU:"PFU",PHIL:"Philips Consumer Electronics",PNTX:"HOYA",POne:"Phase One A/S",PREM:"Premier Computer Innovations",PRIN:"Princeton Graphic Systems",PRIP:"Princeton Publishing Labs",QLUX:"Hong Kong",QMS:"QMS",QPCD:"QPcard AB",QUAD:"QuadLaser",QUME:"Qume",RADI:"Radius",RDDx:"Integrated Color Solutions",RDG:"Roland DG",REDM:"REDMS Group",RELI:"Relisys",RGMS:"Rolf Gierling Multitools",RICO:"Ricoh",RNLD:"Edmund Ronald",ROYA:"Royal",RPC:"Ricoh Printing Systems",RTL:"Royal Information Electronics",SAMP:"Sampo",SAMS:"Samsung",SANT:"Jaime Santana Pomares",SCIT:"Scitex",SCRN:"Dainippon Screen",SDP:"Scitex",SEC:"Samsung",SEIK:"Seiko Instruments",SEIk:"Seikosha",SGUY:"ScanGuy.com",SHAR:"Sharp Laboratories",SICC:"International Color Consortium",SONY:"Sony",SPCL:"SpectraCal",STAR:"Star",STC:"Sampo Technology",Scit:"Scitex",Sdp:"Scitex",Sony:"Sony",TALO:"Talon Technology",TAND:"Tandy",TATU:"Tatung",TAXA:"TAXAN America",TDS:"Tokyo Denshi Sekei",TECO:"TECO Information Systems",TEGR:"Tegra",TEKT:"Tektronix",TI:"Texas Instruments",TMKR:"TypeMaker",TOSB:"Toshiba",TOSH:"Toshiba",TOTK:"TOTOKU ELECTRIC",TRIU:"Triumph",TSBT:"Toshiba",TTX:"TTX Computer Products",TVM:"TVM Professional Monitor",TW:"TW Casper",ULSX:"Ulead Systems",UNIS:"Unisys",UTZF:"Utz Fehlau & Sohn",VARI:"Varityper",VIEW:"Viewsonic",VISL:"Visual communication",VIVO:"Vivo Mobile Communication",WANG:"Wang",WLBR:"Wilbur Imaging",WTG2:"Ware To Go",WYSE:"WYSE Technology",XERX:"Xerox",XRIT:"X-Rite",ZRAN:"Zoran",Zebr:"Zebra Technologies",appl:"Apple Computer",bICC:"basICColor",berg:"bergdesign",ceyd:"Integrated Color Solutions",clsp:"MacDermid ColorSpan",ds:"Dainippon Screen",dupn:"DuPont",ffei:"FujiFilm Electronic Imaging",flux:"FluxData",iris:"Scitex",kart:"Scitex",lcms:"Little CMS",lino:"Linotronic",none:"none",ob4d:"Erdt Systems",obic:"Medigraph",quby:"Qubyx Sarl",scit:"Scitex",scrn:"Dainippon Screen",sdp:"Scitex",siwi:"SIWI GRAFIKA",yxym:"YxyMaster"},Ir={scnr:"Scanner",mntr:"Monitor",prtr:"Printer",link:"Device Link",abst:"Abstract",spac:"Color Space Conversion Profile",nmcl:"Named Color",cenc:"ColorEncodingSpace profile",mid:"MultiplexIdentification profile",mlnk:"MultiplexLink profile",mvis:"MultiplexVisualization profile",nkpf:"Nikon Input Device Profile (NON-STANDARD!)"};Pt(oe,"icc",[[4,Xn],[12,Ir],[40,Object.assign({},Xn,Ir)],[48,Xn],[80,Xn],[64,{0:"Perceptual",1:"Relative Colorimetric",2:"Saturation",3:"Absolute Colorimetric"}],["tech",{amd:"Active Matrix Display",crt:"Cathode Ray Tube Display",kpcd:"Photo CD",pmd:"Passive Matrix Display",dcam:"Digital Camera",dcpj:"Digital Cinema Projector",dmpc:"Digital Motion Picture Camera",dsub:"Dye Sublimation Printer",epho:"Electrophotographic Printer",esta:"Electrostatic Printer",flex:"Flexography",fprn:"Film Writer",fscn:"Film Scanner",grav:"Gravure",ijet:"Ink Jet Printer",imgs:"Photo Image Setter",mpfr:"Motion Picture Film Recorder",mpfs:"Motion Picture Film Scanner",offs:"Offset Lithography",pjtv:"Projection Television",rpho:"Photographic Paper Printer",rscn:"Reflective Scanner",silk:"Silkscreen",twax:"Thermal Wax Printer",vidc:"Video Camera",vidm:"Video Monitor"}]]);class Gn extends ae{static canHandle(t,e,i){return t.getUint8(e+1)===237&&t.getString(e+4,9)==="Photoshop"&&this.containsIptc8bim(t,e,i)!==void 0}static headerLength(t,e,i){let r,s=this.containsIptc8bim(t,e,i);if(s!==void 0)return r=t.getUint8(e+s+7),r%2!=0&&(r+=1),r===0&&(r=4),s+8+r}static containsIptc8bim(t,e,i){for(let r=0;r<i;r++)if(this.isIptcSegmentHead(t,e+r))return r}static isIptcSegmentHead(t,e){return t.getUint8(e)===56&&t.getUint32(e)===943868237&&t.getUint16(e+4)===1028}parse(){let{raw:t}=this,e=this.chunk.byteLength-1,i=!1;for(let r=0;r<e;r++)if(this.chunk.getUint8(r)===28&&this.chunk.getUint8(r+1)===2){i=!0;let s=this.chunk.getUint16(r+3),o=this.chunk.getUint8(r+2),a=this.chunk.getLatin1String(r+5,s);t.set(o,this.pluralizeValue(t.get(o),a)),r+=4+s}else if(i)break;return this.translate(),this.output}pluralizeValue(t,e){return t!==void 0?t instanceof Array?(t.push(e),t):[t,e]:e}}nt(Gn,"type","iptc"),nt(Gn,"translateValues",!1),nt(Gn,"reviveValues",!1),Ct.set("iptc",Gn),Pt(Nt,"iptc",[[0,"ApplicationRecordVersion"],[3,"ObjectTypeReference"],[4,"ObjectAttributeReference"],[5,"ObjectName"],[7,"EditStatus"],[8,"EditorialUpdate"],[10,"Urgency"],[12,"SubjectReference"],[15,"Category"],[20,"SupplementalCategories"],[22,"FixtureIdentifier"],[25,"Keywords"],[26,"ContentLocationCode"],[27,"ContentLocationName"],[30,"ReleaseDate"],[35,"ReleaseTime"],[37,"ExpirationDate"],[38,"ExpirationTime"],[40,"SpecialInstructions"],[42,"ActionAdvised"],[45,"ReferenceService"],[47,"ReferenceDate"],[50,"ReferenceNumber"],[55,"DateCreated"],[60,"TimeCreated"],[62,"DigitalCreationDate"],[63,"DigitalCreationTime"],[65,"OriginatingProgram"],[70,"ProgramVersion"],[75,"ObjectCycle"],[80,"Byline"],[85,"BylineTitle"],[90,"City"],[92,"Sublocation"],[95,"State"],[100,"CountryCode"],[101,"Country"],[103,"OriginalTransmissionReference"],[105,"Headline"],[110,"Credit"],[115,"Source"],[116,"CopyrightNotice"],[118,"Contact"],[120,"Caption"],[121,"LocalCaption"],[122,"Writer"],[125,"RasterizedCaption"],[130,"ImageType"],[131,"ImageOrientation"],[135,"LanguageIdentifier"],[150,"AudioType"],[151,"AudioSamplingRate"],[152,"AudioSamplingResolution"],[153,"AudioDuration"],[154,"AudioOutcue"],[184,"JobID"],[185,"MasterDocumentID"],[186,"ShortDocumentID"],[187,"UniqueDocumentID"],[188,"OwnerID"],[200,"ObjectPreviewFileFormat"],[201,"ObjectPreviewFileVersion"],[202,"ObjectPreviewData"],[221,"Prefs"],[225,"ClassifyState"],[228,"SimilarityIndex"],[230,"DocumentNotes"],[231,"DocumentHistory"],[232,"ExifCameraInfo"],[255,"CatalogSets"]]),Pt(oe,"iptc",[[10,{0:"0 (reserved)",1:"1 (most urgent)",2:"2",3:"3",4:"4",5:"5 (normal urgency)",6:"6",7:"7",8:"8 (least urgent)",9:"9 (user-defined priority)"}],[75,{a:"Morning",b:"Both Morning and Evening",p:"Evening"}],[131,{L:"Landscape",P:"Portrait",S:"Square"}]]);let Ni=null;async function xa(){if(!Ni)try{const i=(await import("./joraw-1Lq5hXK7.js")).default;if(typeof i!="function")throw new Error("JoRaw WASM import failed");const r=new URL("/assets/joraw-DraTMNgX.wasm",import.meta.url).href;Ni=i({locateFile:(s,o)=>s.endsWith("joraw.wasm")?r:o+s})}catch(e){throw console.error("Failed to load joraw.js:",e),e}const n=await Ni,t=n.LibRaw||n.JoRaw;if(!t)throw new Error("JoRaw class not found");return t}const ba=async n=>{var i,r,s,o,a,l,u,d,c;const t=await xa(),e=new t;try{if(await e.open(n,{}),typeof e.getRawImage!="function")throw new Error("WASM mismatch");const g=e.getRawImage();let h=new Uint16Array(g.data);const f=await e.metadata(!0);let p={...f};try{const S=await pa.parse(n.buffer);S&&(p={...p,...S})}catch(S){console.warn("exifr parsing failed for RAW buffer",S)}const m=((i=f.idata)==null?void 0:i.filters)||0,y=((r=f.idata)==null?void 0:r.colors)||0,x=m===0&&y===3,b=m===9;let M=[0,0,0,0],C=!1;if(e.getBlackLevels)try{const S=e.getBlackLevels();S.dng_cblack&&S.dng_cblack.length===4&&Array.from(S.dng_cblack).some(w=>w>0)?(M=Array.from(S.dng_cblack).map(Number),C=!0):S.cblack&&S.cblack.length===4&&Array.from(S.cblack).some(w=>w>0)?(M=Array.from(S.cblack).map(Number),C=!0):typeof S.black=="number"&&S.black>0&&(M=[S.black,S.black,S.black,S.black],C=!0)}catch(S){console.warn("getBlackLevels binding failed",S)}if(!C){let S=[];if((s=f.color_data)!=null&&s.cblack_rawpy_style)S=f.color_data.cblack_rawpy_style;else if((a=(o=f.color_data)==null?void 0:o.dng_levels)!=null&&a.dng_cblack)S=f.color_data.dng_levels.dng_cblack;else if(((l=f.black_level_per_channel)==null?void 0:l.length)>=4)S=f.black_level_per_channel;else if(((u=f.cblack)==null?void 0:u.length)>=4)S=f.cblack;else if(((c=(d=f.color)==null?void 0:d.cblack)==null?void 0:c.length)>=4)S=f.color.cblack;else{const w=f.black_level||f.color_data&&f.color_data.black||0;S=[w,w,w,w]}M=[Number(S[0])||0,Number(S[1])||0,Number(S[2])||0,Number(S[3])||0]}return{data:h,width:g.width,height:g.height,bayerPattern:f.color_desc||"RGGB",blackLevels:M,whiteLevel:f.white_level||16383,metadata:p,isThreePlane:x,threePlaneTransfer:x?"linear":void 0,isXTrans:b}}finally{e.delete?e.delete():e.close()}};async function Ma(n){if(Vs(n)){const i=await Ws(n);if(!i)throw new Error("Sony cRAW HQ decoder did not return image data.");return i.rawImageData}const e=new Uint8Array(n);return ba(e)}function _a(n,t,e){const i=rs(n,e),r=Math.floor(t.x),s=Math.floor(t.y),o=Math.floor(t.w),a=Math.floor(t.h),l=new Uint16Array(o*a);for(let u=0;u<a;u++){const d=s+u,c=u*o;for(let g=0;g<o;g++)l[c+g]=i(r+g,d)}return{data:l,width:o,height:a}}function wa(n,t,e){const i=Sa(n,t,e);if(!i)return null;const r=n.width,s=n.height,o=new Uint16Array(r*s);for(let a=0;a<s;a++){const l=a*r;for(let u=0;u<r;u++)o[l+u]=i(u,a)}return{kind:"u16-mono",data:o,width:r,height:s}}function Sa(n,t,e){return t.renderMode==="advanced-zero-dep"&&t.advancedZeroDep?rs(n,t,e):t.renderMode==="zero-dependency"?va(n,t,e):null}function rs(n,t,e){if(!t.advancedZeroDep)throw new Error("Unmixing settings not found in DisplaySettings.");const{bg:i,fg:r}=t.advancedZeroDep,s=ss(e,t.advancedZeroDep.bl),{data:o,width:a,whiteLevel:l}=n,u=i.map((f,p)=>Math.max(0,f-s[p])),d=r.map((f,p)=>Math.max(0,f-s[p])),c=(u[1]+u[3])/2,g=(d[1]+d[3])/2,h=Math.pow(2,t.exposure);return(f,p)=>{if(f<0||p<0||f>=a||p>=n.height)return 0;const m=p%2,y=f%2;let x=0;!m&&!y?x=0:!m&&y?x=1:m&&!y?x=3:x=2;const b=o[p*a+f],M=s[x],C=u[x],S=d[x],w=Math.max(b-M,0),P=S-C||1e-9,T=(w-C)/P;let F;return T<0?F=w*(c/Math.max(C,1e-9)):T>1?F=w*(g/Math.max(S,1e-9)):F=(1-T)*c+T*g,F*=h,Math.max(0,Math.min(65535,Math.round(F)))}}function va(n,t,e){const{data:i,width:r,height:s,whiteLevel:o}=n,a=ss(e,t.blackLevel||[0,0,0,0]),l=Ca(n.bayerPattern),u=t.wbGains?t.wbGains[0]:1,d=t.wbGains?t.wbGains[1]:1,c=Math.pow(2,t.exposure||0);return(g,h)=>{if(g<0||h<0||g>=r||h>=s)return 0;const f=as(g,h),p=ka(l,g,h),m=i[h*r+g],y=a[f];let x=(m-y)/Math.max(1,o-y);return x=Math.max(0,Math.min(1,x)),x*=c,p==="R"?x*=u:p==="B"&&(x*=d),Pa(x)}}function ss(n,t){if(typeof n=="number"&&Number.isFinite(n)){const e=Math.max(0,n);return[e,e,e,e]}return Array.isArray(n)&&n.length===4?[Math.max(0,n[0]??0),Math.max(0,n[1]??0),Math.max(0,n[2]??0),Math.max(0,n[3]??0)]:[Math.max(0,t[0]??0),Math.max(0,t[1]??0),Math.max(0,t[2]??0),Math.max(0,t[3]??0)]}function Pa(n){return Math.max(0,Math.min(65535,Math.round(Math.max(0,Math.min(1,n))*65535)))}function Ca(n){const t=(n||"RGGB").toUpperCase().trim();return t.length>=4&&/^[RGB]{4}$/.test(t.slice(0,4))?t.slice(0,4):"RGGB"}function as(n,t){return(t&1)<<1|n&1}function ka(n,t,e){return n[as(t,e)]}const ht=(n,t=0)=>({real:n,imag:t}),hn=(n,t)=>({real:n.real+t.real,imag:n.imag+t.imag}),cn=(n,t)=>({real:n.real-t.real,imag:n.imag-t.imag}),It=(n,t)=>({real:n.real*t.real-n.imag*t.imag,imag:n.real*t.imag+n.imag*t.real}),Ee=(n,t)=>{const e=t.real*t.real+t.imag*t.imag;return e===0?ht(0):{real:(n.real*t.real+n.imag*t.imag)/e,imag:(n.imag*t.real-n.real*t.imag)/e}},Qe=n=>Math.hypot(n.real,n.imag),os=n=>{const t=Qe(n);if(t===0)return ht(0);const e=Math.sqrt(t),i=Math.atan2(n.imag,n.real);return ht(e*Math.cos(i/2),e*Math.sin(i/2))};function Fa(n,t){const e=n.length-1;if(e<0)return{p:ht(0),dp:ht(0),d2p:ht(0)};let i=ht(n[e].real,n[e].imag),r=ht(0),s=ht(0);for(let o=e-1;o>=0;o--)s=hn(It(r,ht(2)),It(t,s)),r=hn(i,It(t,r)),i=hn(ht(n[o].real,n[o].imag),It(t,i));return{p:i,dp:r,d2p:s}}function Xi(n,t,e=80){const r=n.length-1;if(r<=0)return{root:t,iterations:0};if(r===1)return{root:Ee(It(n[0],ht(-1)),n[1]),iterations:0};let s=ht(t.real,t.imag);for(let o=0;o<e;o++){const{p:a,dp:l,d2p:u}=Fa(n,s);if(Qe(a)<1e-14)return{root:s,iterations:o};const d=Ee(l,a),c=It(d,d),g=cn(c,Ee(u,a)),h=ht(r),f=ht(r-1),p=cn(It(h,g),It(d,d)),m=os(It(f,p)),y=hn(d,m),x=cn(d,m),b=Qe(y)>Qe(x)?y:x;if(Qe(b)<1e-14)return{root:s,iterations:o};const M=Ee(h,b),C=cn(s,M);if(Qe(M)<1e-14*Qe(C))return{root:C,iterations:o+1};s=C}return{root:s,iterations:e}}function Ta(n,t){const e=n.length-1;if(e<=0)return[ht(0)];if(e===1)return[ht(n[0].real,n[0].imag)];const i=new Array(e);i[e-1]=ht(n[e].real,n[e].imag);for(let r=e-2;r>=0;r--){const s=ht(n[r+1].real,n[r+1].imag),o=i[r+1];i[r]=hn(s,It(t,o))}return i}function Aa(n){const t=n.length-1;if(t<=0)return[];if(t===1)return[Ee(It(n[0],ht(-1)),n[1])];const e=[];let i=n.map(s=>ht(s.real,s.imag)),r=t*5;for(;i.length>2&&r-- >0;){const s=ht(.3+Math.random()*.7,.3+Math.random()*.7),{root:o}=Xi(i,s,100),a=Xi(n,o,20);e.push(a.root);const l=Ta(i,o);if(l.length>=i.length){console.warn("polyDeflate did not reduce degree, breaking");break}i=l}if(i.length===2)e.push(Ee(It(i[0],ht(-1)),i[1]));else if(i.length===3){const s=i[2],o=i[1],a=i[0],l=cn(It(o,o),It(It(ht(4),s),a)),u=os(l),d=It(ht(2),s),c=Ee(cn(It(ht(-1),o),u),d),g=Ee(hn(It(ht(-1),o),u),d);e.push(c,g)}return e}function Ia(n,t,e){const i=[ht(n),ht(-1),ht(n*t),ht(0),ht(n*e)],r=Aa(i);if(r.length===0)return console.warn("laguerreSmallestPositiveRoot: no roots found"),n;let s=1/0,o=!1;for(const l of r)Math.abs(l.imag)<1e-10&&l.real>0&&l.real<s&&(s=l.real,o=!0);return o?Xi(i,ht(s,0),20).root.real:(console.warn("laguerreSmallestPositiveRoot: no positive real root found"),n)}function Na(n,t,e){if(Math.abs(t)<1e-10&&Math.abs(e)<1e-10)return n;if(n<1e-10)return 0;if(Math.abs(e)<1e-10){const i=-1/(t*n),r=1/t,s=i*i-4*r;if(s<0)return n;const o=Math.sqrt(s),a=-.5*(i+Math.sign(i)*o),l=a,u=r/a;return l>0&&u>0?Math.min(l,u):l>0?l:u>0?u:n}try{return Ia(n,t,e)}catch(i){return console.error("laguerreSmallestPositiveRoot failed:",i),n}}function Ra(n,t,e,i,r){const s=n.x-t.x,o=n.y-t.y,a=Math.hypot(s,o)/Math.max(1e-12,e);if(a<1e-12)return{x:n.x,y:n.y};const l=a*a,u=1+(i+r*l)*l;return{x:s/u+t.x,y:o/u+t.y}}function Ea(n,t,e,i,r){const s=n.x-t.x,o=n.y-t.y,a=Math.hypot(s,o)/Math.max(1e-12,e);if(a<1e-12)return{x:t.x,y:t.y};const u=Na(a,i,r)/a;return{x:t.x+s*u,y:t.y+o*u}}function Dt(n,t){const e=Ra(n,{x:t.principalX,y:t.principalY},t.radiusNorm,t.k1,t.k2);return{x:e.x+(t.correctedOffsetX??0),y:e.y+(t.correctedOffsetY??0)}}function Ce(n,t){const e={x:n.x-(t.correctedOffsetX??0),y:n.y-(t.correctedOffsetY??0)};return Ea(e,{x:t.principalX,y:t.principalY},t.radiusNorm,t.k1,t.k2)}const La=`
attribute vec2 a_position;

void main() {
  gl_Position = vec4(a_position, 0.0, 1.0);
}
`,Da=`
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
`,Ba=`
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
`;class Ua{constructor(){bt(this,"canvas",null);bt(this,"gl",null);bt(this,"blurProgram",null);bt(this,"sobelProgram",null);bt(this,"positionBuffer",null);bt(this,"blurUniforms",null);bt(this,"sobelUniforms",null);bt(this,"resources",null);bt(this,"initialized",!1);bt(this,"unavailable",!1);bt(this,"maxTextureSize",0)}compute(t,e,i){if(!this.initialized&&!this.init())return null;const r=this.gl,s=this.blurProgram,o=this.sobelProgram,a=this.blurUniforms,l=this.sobelUniforms;if(!r||!s||!o||!a||!l||!this.positionBuffer||!this.canvas||e<=2||i<=2||e>this.maxTextureSize||i>this.maxTextureSize)return null;const u=this.ensureResources(e,i);if(!u)return null;this.canvas.width=e,this.canvas.height=i,r.viewport(0,0,e,i),r.disable(r.BLEND),r.pixelStorei(r.UNPACK_ALIGNMENT,1),r.pixelStorei(r.PACK_ALIGNMENT,1),r.bindBuffer(r.ARRAY_BUFFER,this.positionBuffer),r.activeTexture(r.TEXTURE0),r.bindTexture(r.TEXTURE_2D,u.sourceTexture),r.texImage2D(r.TEXTURE_2D,0,r.LUMINANCE,e,i,0,r.LUMINANCE,r.UNSIGNED_BYTE,t),r.useProgram(s),r.enableVertexAttribArray(0),r.vertexAttribPointer(0,2,r.FLOAT,!1,0,0),r.uniform2f(a.size,e,i),r.uniform1i(a.source,0),r.bindFramebuffer(r.FRAMEBUFFER,u.blurFramebuffer),r.bindTexture(r.TEXTURE_2D,u.sourceTexture),r.drawArrays(r.TRIANGLES,0,6);const d=new Uint8Array(e*i*4);r.readPixels(0,0,e,i,r.RGBA,r.UNSIGNED_BYTE,d),r.useProgram(o),r.uniform2f(l.size,e,i),r.uniform1i(l.blurred,0),r.bindFramebuffer(r.FRAMEBUFFER,u.sobelFramebuffer),r.bindTexture(r.TEXTURE_2D,u.blurTexture),r.drawArrays(r.TRIANGLES,0,6);const c=new Uint8Array(e*i*4);r.readPixels(0,0,e,i,r.RGBA,r.UNSIGNED_BYTE,c),r.disableVertexAttribArray(0),r.bindFramebuffer(r.FRAMEBUFFER,null),r.bindBuffer(r.ARRAY_BUFFER,null),r.bindTexture(r.TEXTURE_2D,null);const g=new Uint8Array(e*i);for(let m=0,y=0;m<g.length;m++,y+=4)g[m]=d[y];const h=new Float32Array(e*i),f=new Float32Array(e*i),p=new Float32Array(e*i);for(let m=0,y=0;m<h.length;m++,y+=4){const x=(c[y]|c[y+1]<<8)-32768,b=(c[y+2]|c[y+3]<<8)-32768;h[m]=x,f[m]=b,p[m]=Math.sqrt(x*x+b*b)}return{blurredGray:g,gx:h,gy:f,magnitude:p}}init(){if(this.initialized&&this.gl&&this.blurProgram&&this.sobelProgram)return!0;if(this.unavailable)return!1;const t=this.createCanvas();if(!t)return this.unavailable=!0,!1;const e=t.getContext("webgl",{alpha:!1,antialias:!1,depth:!1,stencil:!1,premultipliedAlpha:!1,preserveDrawingBuffer:!1});if(!e)return this.unavailable=!0,!1;const i=this.compileShader(e,e.VERTEX_SHADER,La),r=this.compileShader(e,e.FRAGMENT_SHADER,Da),s=this.compileShader(e,e.FRAGMENT_SHADER,Ba);if(!i||!r||!s)return i&&e.deleteShader(i),r&&e.deleteShader(r),s&&e.deleteShader(s),this.unavailable=!0,!1;const o=this.createProgram(e,i,r),a=this.createProgram(e,i,s);if(e.deleteShader(i),e.deleteShader(r),e.deleteShader(s),!o||!a)return o&&e.deleteProgram(o),a&&e.deleteProgram(a),this.unavailable=!0,!1;const l=e.createBuffer();return l?(e.bindBuffer(e.ARRAY_BUFFER,l),e.bufferData(e.ARRAY_BUFFER,new Float32Array([-1,-1,1,-1,-1,1,-1,1,1,-1,1,1]),e.STATIC_DRAW),e.bindBuffer(e.ARRAY_BUFFER,null),this.canvas=t,this.gl=e,this.blurProgram=o,this.sobelProgram=a,this.positionBuffer=l,this.blurUniforms={source:e.getUniformLocation(o,"u_source"),size:e.getUniformLocation(o,"u_size")},this.sobelUniforms={blurred:e.getUniformLocation(a,"u_blurred"),size:e.getUniformLocation(a,"u_size")},this.maxTextureSize=Number(e.getParameter(e.MAX_TEXTURE_SIZE)||0),this.initialized=!0,!0):(e.deleteProgram(o),e.deleteProgram(a),this.unavailable=!0,!1)}createCanvas(){return typeof OffscreenCanvas<"u"?new OffscreenCanvas(1,1):typeof document<"u"?document.createElement("canvas"):null}ensureResources(t,e){const i=this.gl;if(!i)return null;if(this.resources&&this.resources.width===t&&this.resources.height===e)return this.resources;this.disposeResources();const r=this.createTexture(i.LUMINANCE,t,e,i.LUMINANCE,i.UNSIGNED_BYTE,null),s=this.createTexture(i.RGBA,t,e,i.RGBA,i.UNSIGNED_BYTE,null),o=this.createTexture(i.RGBA,t,e,i.RGBA,i.UNSIGNED_BYTE,null),a=this.createFramebuffer(s),l=this.createFramebuffer(o);return!r||!s||!o||!a||!l?(r&&i.deleteTexture(r),s&&i.deleteTexture(s),o&&i.deleteTexture(o),a&&i.deleteFramebuffer(a),l&&i.deleteFramebuffer(l),null):(this.resources={width:t,height:e,sourceTexture:r,blurTexture:s,sobelTexture:o,blurFramebuffer:a,sobelFramebuffer:l},this.resources)}createTexture(t,e,i,r,s,o){const a=this.gl;if(!a)return null;const l=a.createTexture();return l?(a.bindTexture(a.TEXTURE_2D,l),a.texParameteri(a.TEXTURE_2D,a.TEXTURE_MIN_FILTER,a.NEAREST),a.texParameteri(a.TEXTURE_2D,a.TEXTURE_MAG_FILTER,a.NEAREST),a.texParameteri(a.TEXTURE_2D,a.TEXTURE_WRAP_S,a.CLAMP_TO_EDGE),a.texParameteri(a.TEXTURE_2D,a.TEXTURE_WRAP_T,a.CLAMP_TO_EDGE),a.texImage2D(a.TEXTURE_2D,0,t,e,i,0,r,s,o),a.bindTexture(a.TEXTURE_2D,null),l):null}createFramebuffer(t){const e=this.gl;if(!e||!t)return null;const i=e.createFramebuffer();if(!i)return null;e.bindFramebuffer(e.FRAMEBUFFER,i),e.framebufferTexture2D(e.FRAMEBUFFER,e.COLOR_ATTACHMENT0,e.TEXTURE_2D,t,0);const r=e.checkFramebufferStatus(e.FRAMEBUFFER);return e.bindFramebuffer(e.FRAMEBUFFER,null),r!==e.FRAMEBUFFER_COMPLETE?(e.deleteFramebuffer(i),null):i}compileShader(t,e,i){const r=t.createShader(e);return r?(t.shaderSource(r,i),t.compileShader(r),t.getShaderParameter(r,t.COMPILE_STATUS)?r:(console.error("[SFR Auto Detect WebGL] shader compile failed",t.getShaderInfoLog(r)),t.deleteShader(r),null)):null}createProgram(t,e,i){const r=t.createProgram();return r?(t.attachShader(r,e),t.attachShader(r,i),t.bindAttribLocation(r,0,"a_position"),t.linkProgram(r),t.getProgramParameter(r,t.LINK_STATUS)?r:(console.error("[SFR Auto Detect WebGL] program link failed",t.getProgramInfoLog(r)),t.deleteProgram(r),null)):null}disposeResources(){const t=this.gl,e=this.resources;if(!t||!e){this.resources=null;return}t.deleteTexture(e.sourceTexture),t.deleteTexture(e.blurTexture),t.deleteTexture(e.sobelTexture),t.deleteFramebuffer(e.blurFramebuffer),t.deleteFramebuffer(e.sobelFramebuffer),this.resources=null}}const Oa=new Ua,Gi={gradientPercentiles:[.82,.88,.92,.95,.98,.995],downsampleMaxSide:1600,minComponentAreaRatio:15e-6,maxComponentAreaRatio:.35,minComponentAreaPx:20,minEdgePoints:24,extentQuantileLow:.02,extentQuantileHigh:.98,cornerTrimRatio:.18,minSpanPx:8,maxAspectRatio:2,bandScale:.16,bandMinPx:1.75,bandMaxPx:14,minPointContrast:6,minSidePoints:3,minCoverageRatio:.15,minCenterCoverageRatio:.2,filterBlockPurity:!0,innerPurityStdScale:1.5,outerMeanSpreadLimit:51,minAxisDot:.6,residualLimitFloor:.01,residualLimitScale:.25,minQuadArea:48,minSideLength:10,minOuterContrast:5,sampleHalfWidthRatio:.25};function Xa(n,t,e,i,r,s){const o=n.width,a=n.height,l=(n.bayerPattern||"RGGB").toUpperCase(),u=[],d=[],c=s!=null&&s.correctedRect?Et*2:Et,g=Math.max(1,Math.min(r,c)),h=e.p2.x-e.p1.x,f=e.p2.y-e.p1.y,p=Math.hypot(h,f);if(!Number.isFinite(p)||p<=1e-6)return null;const m=h/p,y=f/p,x=-y,b=m,M={x:(e.p1.x+e.p2.x)*.5,y:(e.p1.y+e.p2.y)*.5},C=s!=null&&s.correctedRect?_t(s.correctedRect,o,a):_t(Xt(xe(e,c*4+2)??[e.p1,e.p2],2),o,a);if(!C)return null;const S=(s==null?void 0:s.correctedScanlinesOverride)??(s!=null&&s.distortedRect?nr(_t(s.distortedRect,o,a)??s.distortedRect,t,o,a):Ss(C,e,Math.max(1,i),g*4+.5,o,a));if(!S||S.size===0)return null;const w=vs(S,t,o,a);if(w.size===0)return null;const P=!ir(t);for(const[F,k]of w){if(F<0||F>=a)continue;const _=F*o;for(let v=k.start;v<=k.end;v++){if(v<0||v>=o||!vt(v,F,l,s==null?void 0:s.greenPhase))continue;const I={x:v,y:F},E=Dt(I,t);if(!Number.isFinite(E.x)||!Number.isFinite(E.y)||Math.round(E.x)<0||Math.round(E.x)>=o||Math.round(E.y)<0||Math.round(E.y)>=a)continue;const L=E.x-M.x,z=E.y-M.y,U=L*m+z*y;let R=L*x+z*b;if(P){const O=Ps(U,m,y,M,I,t);if(!O)continue;const B=.5*(O.a+O.b),Y=qt(O.a,m,y,M,t),Z=qt(B,m,y,M,t),X=qt(O.b,m,y,M,t),J=sr({x:O.a,y:Math.hypot(Y.x-I.x,Y.y-I.y)},{x:B,y:Math.hypot(Z.x-I.x,Z.y-I.y)},{x:O.b,y:Math.hypot(X.x-I.x,X.y-I.y)});if(!Number.isFinite(J))continue;const W=rr(J,m,y,M,t),tt=Math.hypot(W.x,W.y);if(!Number.isFinite(tt)||tt<=1e-9)continue;const lt=W.x/tt,$=-(W.y/tt),K=lt,A=qt(J,m,y,M,t);R=(I.x-A.x)*$+(I.y-A.y)*K}!Number.isFinite(U)||Math.abs(U)>Math.max(1,i)||!Number.isFinite(R)||Math.abs(R)>g||(u.push(R),d.push(Math.max(0,n.data[_+v]-gn(s==null?void 0:s.blackLevel,v,F))))}}if(u.length<8)return null;const T=Math.abs(h)>=Math.abs(f)?1:2;return s!=null&&s.forceLegacyModel?kn(u,d,T,c):yn(u,d,T,c)}function $i(n,t){const e=n.length;let i=0,r=0,s=0,o=0;for(let l=0;l<e;l++)i+=n[l],r+=t[l],s+=n[l]*t[l],o+=n[l]*n[l];const a=e*o-i*i;return a===0?{slope:0,intercept:0}:{slope:(e*s-i*r)/a,intercept:(r*o-i*s)/a}}function Ga(n,t){const e=n.length,i=new Array(e).fill(0),r=2*t;for(let s=0;s<e;s++){const o=s>0?n[s-1]:n[0],a=s<e-1?n[s+1]:n[e-1];i[s]=(a-o)/r}return i}const Zt=-1e7,dn=13,se=512,Lt=8,li=1/Lt,Et=28,Ya=[[0,0,0,0,0,-.085714285714286,.342857142857143,.485714285714286,.342857142857143,-.085714285714286,0,0,0,0,0],[0,0,0,0,-.095238095238095,.142857142857143,.285714285714286,.333333333333333,.285714285714286,.142857142857143,-.095238095238095,0,0,0,0],[0,0,0,-.090909090909091,.060606060606061,.168831168831169,.233766233766234,.255411255411255,.233766233766234,.168831168831169,.060606060606061,-.090909090909091,0,0,0],[0,0,-.083916083916084,.020979020979021,.102564102564103,.160839160839161,.195804195804196,.207459207459208,.195804195804196,.160839160839161,.102564102564103,.020979020979021,-.083916083916084,0,0],[0,-.076923076923077,0,.062937062937063,.111888111888112,.146853146853147,.167832167832168,.174825174825175,.167832167832168,.146853146853147,.111888111888112,.062937062937063,0,-.076923076923077,0],[-.070588235294118,-.011764705882353,.038009049773756,.078733031674208,.110407239819004,.133031674208145,.146606334841629,.151131221719457,.146606334841629,.133031674208145,.110407239819004,.078733031674208,.038009049773756,-.011764705882353,-.070588235294118]];function fn(n,t,e,i=1){const r=Math.max(1e-6,e*.5),s=Math.max(1e-6,t*i),o=Math.exp(-s*r),a=1-o;if(!Number.isFinite(a)||Math.abs(a)<=1e-9)return Math.abs(n)<=r?1:0;if(Math.abs(n)<r){const l=2-2*o*Math.cosh(s*n),u=2*Math.sinh(s*r)*a;return!Number.isFinite(u)||Math.abs(u)<=1e-9?0:l/u}return Math.exp(-s*Math.abs(n))/a}function ls(n,t,e,i,r,s){const o=n.length;if(o===0)return[];if(s<1)return n;const a=Math.min(s,32),l=new Array(o).fill(0);l[0]=n[0];for(let h=1;h<o;h++)l[h]=l[h-1]+n[h];const u=(h,f)=>{const p=Math.max(0,h),m=Math.min(o-1,f);return m<p?n[Math.max(0,Math.min(o-1,h))]??0:(l[m]-(p>0?l[p-1]:0))/(m-p+1)},d=a*2,c=a,g=1;for(let h=Math.max(t+c,i-d);h<i;h++){const f=Math.max(g,Math.trunc((i-h)*c/Math.max(1,d)));n[h]=u(h-f,h+f)}for(let h=Math.min(r+d-1,e-c-1);h>r;h--){const f=Math.max(g,Math.trunc((h-r)*c/Math.max(1,d)));n[h]=u(h-f,h+f)}for(let h=c+1;h<i-d;h++)n[h]=u(h-c,h+c);for(let h=Math.min(r+d,e-c-1);h<o-c-1;h++)n[h]=u(h-c,h+c);return n}function cs(n){return!Number.isFinite(n)||Math.abs(n)<=1e-9?1:Math.sin(n)/n}let Yn=null,Ri=null;function za(){if(Yn)return Yn;const n=.625,t=1/128,e=Math.max(16,Math.round(n*2/t)+1),i=[],r=[];for(let c=0;c<e;c++){const g=-n+c*t;i.push(g),r.push(Math.abs(g)<=n?fn(g,dn,.125,1):0)}const s=4,o=1/1024,a=Math.round(s/o)+1,l=new Array(a).fill(0),u=new Array(a).fill(1);let d=0;for(let c=0;c<r.length;c++)d+=r[c];d=Math.max(1e-9,d);for(let c=0;c<a;c++){const g=c*o;l[c]=g;let h=0;for(let f=0;f<i.length;f++)h+=r[f]*Math.cos(2*Math.PI*g*i[f]);u[c]=Math.max(1e-6,Math.abs(h)/d)}return Yn={freqs:l,values:u},Yn}function Va(n,t){const e=Math.max(1e-6,cs(Math.PI*n*t)),i=za(),r=Je(Math.max(0,Math.min(i.freqs[i.freqs.length-1],n)),i.freqs,i.values);return Math.max(1e-6,e*r)}function ja(){if(Ri)return Ri;const n=new Array(se/16*4).fill(1),t=se*16,e=new Float32Array(t);for(let s=0;s<t;s++){const o=(s-t/2)/(16*Lt);e[s]=Math.abs(o)<=.625?fn(o,dn,li,1):0}const i=new Cn(t);i.transform(e);const r=Math.max(1e-9,Math.abs(i._real[0]));n[0]=1;for(let s=1;s<n.length;s++){const o=cs(Math.PI*s/256),a=Math.max(1e-6,Math.hypot(i._real[s],i._imag[s])/r);n[s]=Math.max(1e-6,o*a)}return Ri=n,n}function Wa(n,t,e){if(n.length===0||t.length!==n.length||!(e>0))return null;const i=n[0],r=n[n.length-1],s=Math.floor((r-i)/e);if(s<16)return null;const o=Math.max(0,Math.min(s-1,Math.round(-i/e))),a=Math.max(2,Math.round(2/e)),l=Math.max(1,Math.round(.5/e)),u=5,d=.125/Math.max(e,1e-6),c=D=>i+D*e,h=((D,$)=>{let K=Math.max(0,D),A=Math.min(s-1,$);if(A-K<8)return null;let Q=0;for(;;){const G=new Array(s).fill(0),V=new Array(s).fill(0);for(let q=0;q<n.length;q++){const Ft=Math.max(0,Math.min(s-1,Math.trunc((n[q]-i)/e))),At=Math.max(K,Ft-u),xt=Math.min(A-1,Ft+u);for(let wt=At;wt<=xt;wt++){const ee=c(wt),$t=Math.max(0,1-Math.abs((n[q]-ee)*1.75*d));$t<=0||(G[wt]+=t[q]*$t,V[wt]+=$t)}}const rt=new Array(s).fill(Zt);let st=0,it=0,dt=0,ft=0,et=-1,Bt=-1;const Gt=Math.max(o-Math.round(s/8),K+a),j=Math.min(o+Math.round(s/8),A-a);for(let q=Math.max(0,K-1);q<=Math.min(s-1,A+1);q++)V[q]>0&&(rt[q]=G[q]/V[q],q<Gt&&(st+=rt[q],dt++),q>j&&(it+=rt[q],ft++),et<0&&(et=q),Bt=q);if(et<0||Bt<0||dt===0||ft===0)return null;for(let q=et-1;q>=0;q--)rt[q]=rt[et];for(let q=Bt+1;q<s;q++)rt[q]=rt[Bt];const ct=Math.max(2,a),at=new Array(s).fill(0);let pt=o;for(let q=ct+1;q<s-1-ct;q++){let Ft=0,At=0;for(let xt=-ct;xt<=ct;xt++)At+=rt[q+xt]*xt,Ft+=xt*xt;at[q]=Ft>0?At/Ft:0,Math.abs(at[q])>Math.abs(at[pt]??0)&&q>K+ct&&q<A-ct-1&&(pt=q)}const mt=Math.max(1,Math.round(2/e)),te=Math.max(mt+1,Math.round(12/e)),Ut=Math.abs(pt-o);if(Ut>mt&&Ut<te)return null;let Mt=0;for(let q=Math.max(0,o-ct);q<=Math.min(s-1,o+ct);q++)Math.abs(at[q])>Math.abs(Mt)&&(Mt=at[q]);if(!Number.isFinite(Mt)||Math.abs(Mt)<=1e-9)return null;const le=Math.abs(Mt*.001);let Kt=!1,ke=!1;for(let q=Math.max(0,o-a);q>=K+l;q--)if(at[q]*Mt<0&&Math.abs(at[q])>le){let Ft=0,At=0,xt=0;for(let wt=q;wt>=K;wt--)at[wt]*Mt<0&&(Ft++,At=Math.max(At,Math.abs(at[wt]))),xt++;if(Ft>xt*.4&&At/Math.abs(Mt)>.25||Ft>.9*xt&&xt>a){K=Math.min(q,o-a),Kt=!0;break}}for(let q=Math.min(s-1,o+a);q<A-l;q++)if(at[q]*Mt<0&&Math.abs(at[q])>le){let Ft=0,At=0,xt=0;for(let wt=q;wt<A;wt++)at[wt]*Mt<0&&(Ft++,At=Math.max(At,Math.abs(at[wt]))),xt++;if(Ft>xt*.4&&At/Math.abs(Mt)>.25||Ft>.9*xt&&xt>a){A=Math.max(q,o+a),Kt=!0;break}}if(Kt&&((o-K<Math.max(1,Math.round(4/e))||A-o<Math.max(1,Math.round(4/e)))&&(ke=!0),Q<2)){Q++;continue}return{sampled:rt,fftLeft:K,fftRight:A,leftMean:st/dt,rightMean:it/ft,peakSlopeIdx:pt,slopes:at,clipped:Kt,dodgy:ke}}})(0,s-1);if(!h)return null;const f=2,p=Math.max(h.leftMean,h.rightMean),m=Math.min(h.leftMean,h.rightMean);let y=h.fftLeft,x=h.fftLeft,b=1/0,M=1/0;for(let D=h.fftLeft;D<=h.fftRight;D++){let $=0,K=0;for(let V=-f;V<=f;V++){const rt=h.sampled[Math.max(0,Math.min(s-1,D+V))];$+=rt,K++}const A=$/Math.max(1,K),Q=Math.abs(A-m-.1*(p-m)),G=Math.abs(A-m-.9*(p-m));Q<b&&(b=Q,y=D),G<M&&(M=G,x=D)}const C=Math.max(4,Math.abs(y-x)*e);if(y<x){const D=y;y=x,x=D}const S=Math.max(l,l+2*Math.round(C/Math.max(e,1e-6)));y+=S,x-=S;const w=Math.max(Math.abs(y-o),Math.abs(x-o),Math.max(a,Math.round(4/Math.max(e,1e-6)))),P=1.85,T=.5,F=new Array(s).fill(0),k=new Array(s).fill(0);for(let D=0;D<n.length;D++){const $=Math.max(0,Math.min(s-1,Math.trunc((n[D]-i)/e)));let K=5;Math.abs($-o)>T*w&&(K=Math.abs($-o)>2*T*w?12:7);const A=Math.max(h.fftLeft,$-K),Q=Math.min(h.fftRight-1,$+K);if(Q<o-P*w||A>o+P*w){for(let G=A;G<=Q;G++)F[G]+=t[D],k[G]+=1;continue}for(let G=A;G<=Q;G++){let V=1;if(Math.abs(G-o)<P*w){const rt=c(G);if(Math.abs(G-o)<w*T)V=fn(n[D]-rt,dn,e,1);else{const st=(Math.abs(G-o)/Math.max(1e-6,w)-T)/Math.max(1e-6,P-T),it=1*(1-st)+.01*st;V=fn(n[D]-rt,dn,e,it)}}!(V>0)||!Number.isFinite(V)||(F[G]+=t[D]*V,k[G]+=V)}}const _=new Array(s).fill(0);let v=-1,I=-1;for(let D=0;D<s;D++)k[D]>0?(_[D]=F[D]/k[D],v<0&&(v=D),I=D):_[D]=Zt;if(v<0||I<0)return null;const E=Math.max(1,Math.round(3/Math.max(e,1e-6)));let L=_[v],z=1;for(let D=v+1;D<o&&z<E;D++)_[D]!==Zt&&(L+=_[D],z++);L/=z;let U=_[I],R=1;for(let D=I-1;D>o&&R<E;D--)_[D]!==Zt&&(U+=_[D],R++);U/=R;for(let D=v-1;D>=0;D--)_[D]=L;for(let D=I+1;D<s;D++)_[D]=U;for(let D=v+1;D<I;D++){if(_[D]!==Zt)continue;let $=D-1;for(;$>=0&&_[$]===Zt;)$--;let K=D+1;for(;K<s&&_[K]===Zt;)K++;if($>=0&&K<s){const A=(D-$)/Math.max(1,K-$);_[D]=_[$]*(1-A)+_[K]*A}else $>=0?_[D]=_[$]:K<s&&(_[D]=_[K])}const B=[...L<=U?_:[..._].reverse()],Y=L<=U?B:B.reverse(),Z=Math.max(Math.round(o-P*w),h.fftLeft+2),X=Math.min(Math.round(o+P*w),h.fftRight-3),J=Math.max(1,Math.round(2/Math.max(e,1e-6))),W=ls(Y,h.fftLeft,h.fftRight,Z,X,J),tt=new Array(s).fill(0);let lt=W[Math.max(0,Math.min(s-1,h.fftLeft))]??W[0]??0;for(let D=h.fftLeft;D<=h.fftRight;D++){const $=W[D]??lt,K=W[Math.min(s-1,D+1)]??$;tt[D]=K-lt,lt=$}return{esf:W,lsfFull:tt}}function us(n,t,e=Et){if(n.length===0||t.length!==n.length)return null;const i=se,r=i/2,s=li,o=2*Lt,a=Math.max(1,Math.round(.5*Lt)),l=5,u=Math.max(0,Math.round(r-e*Lt)),d=Math.min(i-1,Math.round(r+e*Lt));if(d-u<32)return null;let c=new Array(i).fill(Zt),g=0,h=0,f=0,p=0,m=u,y=d,x=u,b=d,M=0;for(;;){const A=new Array(i).fill(0),Q=new Array(i).fill(0);c=new Array(i).fill(Zt),g=0,h=0,f=0,p=0;let G=-1,V=-1;for(let j=0;j<n.length;j++){const ct=Math.trunc(n[j]*Lt+r),at=Math.max(x,ct-l),pt=Math.min(b-1,ct+l);for(let mt=at;mt<=pt;mt++){const te=(mt-r)*s,Ut=Math.max(0,1-Math.abs((n[j]-te)*1.75));Ut>0&&(Q[mt]+=t[j]*Ut,A[mt]+=Ut)}}const rt=Math.max(r-i/8,x+2*Lt),st=Math.min(r+i/8,b-2*Lt);for(let j=Math.max(0,x-1);j<=Math.min(i-1,b+1);j++)A[j]>0&&(c[j]=Q[j]/A[j],j<rt&&(g+=c[j],f++),j>st&&(h+=c[j],p++),G<0&&(G=j),V=j);if(G<0||V<0||f<=0||p<=0)return null;for(let j=G-1;j>=0;j--)c[j]=c[G];for(let j=V+1;j<i;j++)c[j]=c[V];const it=new Array(i).fill(0);let dt=r;const ft=2*Lt;for(let j=ft+1;j<i-1-ft;j++){let ct=0,at=0;for(let pt=-ft;pt<=ft;pt++)at+=c[j+pt]*pt,ct+=pt*pt;it[j]=ct>0?at/ct:0,Math.abs(it[j])>Math.abs(it[dt]??0)&&j>x+ft&&j<b-ft-1&&(dt=j)}if(Math.abs(dt-r)>2*Lt&&Math.abs(dt-r)<12*Lt)return null;let et=0;for(let j=Math.max(0,r-ft);j<=Math.min(i-1,r+ft);j++)Math.abs(it[j])>Math.abs(et)&&(et=it[j]);if(!Number.isFinite(et)||Math.abs(et)<=1e-9)return null;const Bt=Math.abs(et*.001);m=x,y=b;let Gt=!1;for(let j=r-o;j>=x+a;j--)if(it[j]*et<0&&Math.abs(it[j])>Bt){let ct=0,at=0,pt=0;for(let mt=j;mt>=x;mt--)it[mt]*et<0&&(ct++,at=Math.max(at,Math.abs(it[mt]))),pt++;if(ct>pt*.4&&at/Math.abs(et)>.25||ct>.9*pt&&pt>o){m=Math.min(j,r-o),Gt=!0;break}}for(let j=r+o;j<b-a;j++)if(it[j]*et<0&&Math.abs(it[j])>Bt){let ct=0,at=0,pt=0;for(let mt=j;mt<b;mt++)it[mt]*et<0&&(ct++,at=Math.max(at,Math.abs(it[mt]))),pt++;if(ct>pt*.4&&at/Math.abs(et)>.25||ct>.9*pt&&pt>o){y=Math.max(j,r+o),Gt=!0;break}}if(Gt&&M<2){x=m,b=y,M++;continue}break}const C=Math.max(g/f,h/p),S=Math.min(g/f,h/p);let w=m,P=m,T=1/0,F=1/0;for(let A=m;A<=y;A++){const Q=(c[Math.max(0,A-2)]+c[Math.max(0,A-1)]+c[A]+c[Math.min(i-1,A+1)]+c[Math.min(i-1,A+2)])/5,G=Math.abs(Q-S-.1*(C-S)),V=Math.abs(Q-S-.9*(C-S));G<T&&(T=G,w=A),V<F&&(F=V,P=A)}if(w<P){const A=w;w=P,P=A}const k=Math.max(4,Math.abs(w-P)*s),_=Math.max(a,a+2*Math.trunc(k/Math.max(s,1e-6)));w+=_,P-=_;const v=Math.max(Math.abs(w-r),Math.abs(P-r),Math.max(o,Math.trunc(4/Math.max(s,1e-6)))),I=new Array(i).fill(0),E=new Array(i).fill(0),L=1.85,z=.5;for(let A=0;A<n.length;A++){const Q=Math.trunc(n[A]*Lt+r);let G=5;Math.abs(Q-r)>z*v&&(G=Math.abs(Q-r)>2*z*v?12:7);const V=Math.max(m,Q-G),rt=Math.min(y-1,Q+G);if(rt<r-L*v||V>r+L*v){for(let st=V;st<=rt;st++)I[st]+=t[A],E[st]+=1;continue}for(let st=V;st<=rt;st++){let it=1;if(Math.abs(st-r)<L*v){const dt=(st-r)*s;if(Math.abs(st-r)<v*z)it=fn(n[A]-dt,dn,s,1);else{const ft=(Math.abs(st-r)/Math.max(1e-6,v)-z)/Math.max(1e-6,L-z),et=1*(1-ft)+.01*ft;it=fn(n[A]-dt,dn,s,et)}}!(it>0)||!Number.isFinite(it)||(I[st]+=t[A]*it,E[st]+=it)}}const U=new Array(i).fill(0);let R=-1,O=-1;for(let A=Math.max(0,m-1);A<=Math.min(i-1,y+1);A++)E[A]>0?(U[A]=I[A]/E[A],R<0&&(R=A),O=A):U[A]=Zt;if(R<0||O<0)return null;const B=3*Lt;let Y=U[R],Z=1;for(let A=R+1;A<r&&Z<B;A++)U[A]!==Zt&&(Y+=U[A],Z++);Y/=Math.max(1,Z);let X=U[O],J=1;for(let A=O-1;A>r&&J<B;A--)U[A]!==Zt&&(X+=U[A],J++);X/=Math.max(1,J);for(let A=R-1;A>=0;A--)U[A]=Y;for(let A=O+1;A<i;A++)U[A]=X;const W=Math.max(Math.trunc(r-L*v),m+2),tt=Math.min(Math.trunc(r+L*v),y-3),lt=Math.max(1,Math.trunc(2/Math.max(s,1e-6))),D=ls(U,m,y,W,tt,lt),$=new Array(i).fill(0);let K=D[Math.max(0,Math.min(i-1,m))]??D[0]??0;for(let A=m;A<=y;A++){const Q=D[A]??K;$[A]=(D[Math.min(i-1,A+1)]??Q)-K,K=Q}return{esf:D,lsfFull:$}}function Ha(n){const t=new Array(n.length).fill(0);if(n.length===0)return t;t[0]=n[0];for(let e=1;e<n.length;e++){const i=Ji(n[e]-n[e-1]);t[e]=t[e-1]+i}return t}function Ji(n){if(!Number.isFinite(n))return n;let t=(n+Math.PI)%(2*Math.PI);return t<0&&(t+=2*Math.PI),t-Math.PI}function Qa(n,t,e=0){if(n.length===0)return[];const i=Number.isFinite(e)?e:0,r=n.map((o,a)=>{const l=t[a]??0,u=-2*Math.PI*i*l;return Ji(o-u)});return Ha(r).map((o,a)=>{const l=t[a]??0,u=-2*Math.PI*i*l;return o+u})}function qa(n,t,e,i=.05,r=Number.POSITIVE_INFINITY){const s=Math.min(n.length,t.length);if(s<2)return null;let o=0;if(e)for(let p=1;p<s;p++){const m=e[p];Number.isFinite(m)&&m>o&&(o=m)}const a=e&&o>0?Math.max(1e-6,o*i):0;let l=0,u=0,d=0,c=0,g=0,h=0;for(let p=1;p<s;p++){const m=n[p],y=t[p];if(!Number.isFinite(m)||!Number.isFinite(y)||Math.abs(m)<=1e-12||m>r)continue;const x=e?e[p]:1;if(!Number.isFinite(x)||x<=a)continue;const b=e?x*x:1;l+=b,u+=b*m,d+=b*y,c+=b*m*m,g+=b*m*y,h++}if(h<4||l<=0)return null;const f=l*c-u*u;return Math.abs(f)<=1e-12?null:{slope:(l*g-u*d)/f,intercept:(d*c-u*g)/f,used:h,threshold:a}}function Ka(n,t,e=Number.POSITIVE_INFINITY){const i=[],r=[],s=Math.min(n.length,t.length);for(let o=1;o<s;o++)Number.isFinite(n[o])&&Number.isFinite(t[o])&&Math.abs(n[o])>1e-12&&n[o]<=e&&(i.push(n[o]),r.push(t[o]));return i.length<2?{slope:0,intercept:Number.isFinite(t[0])?t[0]:0,used:i.length}:{...$i(i,r),used:i.length}}function $a(n,t,e,i=.05,r=Number.POSITIVE_INFINITY,s=0){const o=Math.min(n.length,t.length);if(o<4)return null;let a=0;if(e)for(let b=1;b<o;b++){const M=e[b];Number.isFinite(M)&&M>a&&(a=M)}const l=e&&a>0?Math.max(1e-6,a*i):0,u=[];for(let b=1;b<o;b++){const M=n[b],C=t[b];if(!Number.isFinite(M)||!Number.isFinite(C)||Math.abs(M)<=1e-12||M>r)continue;const S=e?e[b]:1;!Number.isFinite(S)||S<=l||u.push({freq:M,phase:C,weight:e?S*S:1})}if(u.length<4)return null;const d=b=>{let M=0,C=0;for(const F of u){const k=F.phase+2*Math.PI*b*F.freq;M+=F.weight*Math.sin(k),C+=F.weight*Math.cos(k)}const S=Math.atan2(M,C);let w=0,P=0;const T=.65;for(const F of u){const k=F.phase+2*Math.PI*b*F.freq,_=Math.abs(Ji(k-S)),v=_<=T?_*_:T*(2*_-T);w+=F.weight*v,P+=F.weight}return{score:P>0?w/P:Number.POSITIVE_INFINITY,intercept:S}},c=Number.isFinite(s)?s:0,g=Math.max(2,Math.min(8,Math.abs(c)>1e-6?4:2)),h=.02;let f=c,p=d(f);for(let b=c-g;b<=c+g+h*.5;b+=h){const M=d(b);M.score<p.score&&(p=M,f=b)}let m=f-h*2,y=f+h*2;for(let b=0;b<32;b++){const M=m+(y-m)/3,C=y-(y-m)/3,S=d(M).score,w=d(C).score;S<w?y=C:m=M}const x=(m+y)*.5;return p=d(x),{slope:-2*Math.PI*x,intercept:p.intercept,used:u.length,threshold:l}}function Ja(n,t){if(Number.isFinite(t)&&t>0)return t;let e=0;for(const i of n)Number.isFinite(i)&&i>e&&(e=i);return e>0?e:Number.POSITIVE_INFINITY}function hs(n,t,e,i=Number.POSITIVE_INFINITY,r=0){const s=Qa(n,t,r),o=Ja(t,i),a=$a(t,n,e,.05,o,r),l=a?null:qa(t,s,e,.05,o),u=a||l?null:Ka(t,s,o),d=(a==null?void 0:a.slope)??(l==null?void 0:l.slope)??(u==null?void 0:u.slope)??0,c=(a==null?void 0:a.intercept)??(l==null?void 0:l.intercept)??(u==null?void 0:u.intercept)??0,g=s.map((x,b)=>x-(d*(t[b]??0)+c)),h=Number.isFinite(g[0])?g[0]:0,f=t.map(x=>d*x+c+h),p=g.map(x=>x-h),m=c+h,y=Number.isFinite(d)?-d/(2*Math.PI):null;return{raw:[...n],unwrapped:s,linear:f,residual:p,fit:{groupDelayPx:y===null?null:y-r,absoluteGroupDelayPx:y,referenceDelayPx:r,slopeRadPerCycle:Number.isFinite(d)?d:null,interceptRad:Number.isFinite(m)?m:null,fitPointCount:(a==null?void 0:a.used)??(l==null?void 0:l.used)??(u==null?void 0:u.used)??0,fitWeightThreshold:(a==null?void 0:a.threshold)??(l==null?void 0:l.threshold)??0,fitDomain:"cycles-per-pixel",fitMaxFreqCyclesPerPixel:o}}}function ds(n,t,e){const i=[],r=[],s=[],o=[];for(let a=0;a<e.length;a++){const l=e[a];i.push(Je(l,t,n.raw)),r.push(Je(l,t,n.unwrapped)),s.push(Je(l,t,n.linear)),o.push(Je(l,t,n.residual))}return{ptfRaw:i,ptfUnwrapped:r,ptfLinear:s,ptfResidual:o}}function fs(n,t){const e=n.map((a,l)=>({dist:a,value:t[l]})).filter(a=>Number.isFinite(a.dist)&&Number.isFinite(a.value)).sort((a,l)=>a.dist-l.dist);if(e.length===0)return{dists:[],vals:[]};const i=Math.max(1,Math.min(16,Math.floor(e.length*.1)));let r=0,s=0;for(let a=0;a<i;a++)r+=e[a].value,s+=e[e.length-1-a].value;r/=i,s/=i;const o=r<=s?e:e.map(a=>({dist:-a.dist,value:a.value})).sort((a,l)=>a.dist-l.dist);return{dists:o.map(a=>a.dist),vals:o.map(a=>a.value)}}function Za(n,t,e="RGGB"){const i=e.toUpperCase(),r=t%2,s=n%2;return i==="RGGB"||i==="BGGR"?(r+s)%2!==0:i==="GBRG"||i==="GRBG"?(r+s)%2===0:(r+s)%2!==0}function to(n,t){return n+t&1?2:1}function Nr(n,t){return(t&1)<<1|n&1}function gn(n,t,e){return n===void 0?0:typeof n=="number"?Number.isFinite(n)?n:0:Number.isFinite(n[Nr(t,e)])?n[Nr(t,e)]:0}function vt(n,t,e,i){return i!==void 0&&i!=="default"?to(n,t)===i:Za(n,t,e)}function Zn(n){return n.length===0?0:n.reduce((t,e)=>t+e,0)/n.length}function eo(n,t,e){const i=(e%t+t)%t,r=Math.floor(i),s=(r+1)%t,o=i-r,a=r<n.length?n[r]:0,l=s<n.length?n[s]:0;return a*(1-o)+l*o}function ps(n,t){const e=n.length;if(e===0)return 0;if(t<=0)return n[0];if(t>=e-1)return n[e-1];const i=Math.floor(t),r=Math.min(e-1,i+1),s=t-i;return n[i]*(1-s)+n[r]*s}function no(n,t,e){const i=n.length,r=new Array(i).fill(0);for(let s=0;s<i;s++)r[s]=ps(n,s-e+t);return r}function ti(n,t,e){const i=n.length;if(i===0)return{peakPos:0,peakIdx:0,peakVal:0};const r=Math.max(0,Math.floor(t-e)),s=Math.min(i-1,Math.ceil(t+e));let o=Math.max(0,Math.min(i-1,Math.round(t))),a=-1/0;for(let u=r;u<=s;u++){const d=Math.abs(n[u]);d>a&&(a=d,o=u)}Number.isFinite(a)||(a=Math.abs(n[o]??0));let l=o;if(o>0&&o<i-1){const u=n[o]>=0?1:-1,d=u*n[o-1],c=u*n[o],g=u*n[o+1],h=d-2*c+g;if(Number.isFinite(h)&&Math.abs(h)>1e-9){const f=.5*(d-g)/h;Number.isFinite(f)&&Math.abs(f)<=1&&(l=o+f)}}return{peakPos:l,peakIdx:o,peakVal:Math.abs(ps(n,l))}}function io(n,t,e,i,r,s,o,a,l){const u=Math.floor(i.x),d=Math.floor(i.y),c=Math.floor(i.w),g=Math.floor(i.h),h=[],f=(p,m)=>{if(p<0||m<0||p>=t||m>=e)return null;const y=r+p,x=s+m;return Math.max(0,n[m*t+p]-gn(l,y,x))};for(let p=0;p<g;p++){const m=[],y=d+p;for(let x=0;x<c;x++){const b=u+x,M=r+b,C=s+y;if(vt(M,C,o,a)){m.push(f(b,y)??0);continue}const S=[],w=f(b-1,y),P=f(b+1,y),T=f(b,y-1),F=f(b,y+1);if(w!==null&&vt(M-1,C,o,a)&&S.push(w),P!==null&&vt(M+1,C,o,a)&&S.push(P),T!==null&&vt(M,C-1,o,a)&&S.push(T),F!==null&&vt(M,C+1,o,a)&&S.push(F),S.length===0){const k=[],_=f(b-1,y-1),v=f(b+1,y-1),I=f(b-1,y+1),E=f(b+1,y+1);_!==null&&vt(M-1,C-1,o,a)&&k.push(_),v!==null&&vt(M+1,C-1,o,a)&&k.push(v),I!==null&&vt(M-1,C+1,o,a)&&k.push(I),E!==null&&vt(M+1,C+1,o,a)&&k.push(E),m.push(Zn(k));continue}m.push(Zn(S))}h.push(m)}return h}function ms(n,t,e,i,r,s,o,a,l){const u=Math.floor(i.x),d=Math.floor(i.y),c=Math.floor(i.w),g=Math.floor(i.h),h=[];for(let f=0;f<g;f++){const p=d+f;for(let m=0;m<c;m++){const y=u+m,x=r+y,b=s+p;vt(x,b,o,a)&&h.push({x,y:b,value:Math.max(0,n[p*t+y]-gn(l,x,b))})}}return h}function gs(n,t,e,i,r){const s=Math.floor(i.x),o=Math.floor(i.y),a=Math.floor(i.w),l=Math.floor(i.h),u=(r==null?void 0:r.globalX)??0,d=(r==null?void 0:r.globalY)??0,c=!!(r!=null&&r.isThreePlane)&&n.length>=t*e*3,g=r==null?void 0:r.threePlaneChannel,h=[];for(let f=0;f<l;f++){const p=o+f,m=p*t;for(let y=0;y<a;y++){const x=s+y;let b=0;if(!c)b=Math.max(0,n[m+x]-gn(r==null?void 0:r.blackLevel,u+x,d+p));else{const M=(m+x)*3;if(g!==void 0)b=n[M+g];else{const C=n[M],S=n[M+1],w=n[M+2];b=.2126*C+.7152*S+.0722*w}}h.push({x:u+x,y:d+p,value:b})}}return h}function ro(n){var s;const t=n.length,e=((s=n[0])==null?void 0:s.length)??0;let i=0,r=0;for(let o=1;o<t-1;o++)for(let a=1;a<e-1;a++)i+=Math.abs(n[o][a+1]-n[o][a-1]),r+=Math.abs(n[o+1][a]-n[o-1][a]);return{gx:i,gy:r}}function ys(n,t,e,i,r,s,o){var g;const a=n.length,l=((g=n[0])==null?void 0:g.length)??0,u=(h,f,p,m,y,x,b)=>{const M=b?a:l,C=Math.max(0,f-3),S=Math.min(M,f+4);let w=0,P=0;for(let F=C;F<S;F++)w+=h[F],P+=F*h[F];if(w<=0)return null;const T=P/w;return b?{x:t+m*y,y:p+T*x,weight:w}:{x:p+T*x,y:e+m*y,weight:w}},d=(h,f,p,m,y,x,b)=>{const M=Math.max(3,Math.min(Math.max(3,Math.floor(p/3)),Math.max(4,Math.round(p*.12)))),C=h.map(v=>{let I=-1/0,E=-1;for(let L=0;L<v.length;L++)v[L]>I&&(I=v[L],E=L);return{peakValue:I,peakIndex:E}}),S=(f-1)*.5,w=C.map((v,I)=>({...v,index:I})).filter(v=>v.peakValue>1&&v.peakIndex>=0).sort((v,I)=>{const E=I.peakValue-v.peakValue;return Math.abs(E)>1e-6?E:Math.abs(v.index-S)-Math.abs(I.index-S)});if(w.length===0)return[];const P=w[0],T=new Array(f).fill(null),F=u(h[P.index],P.peakIndex,m,P.index,y,x,b);if(!F)return[];T[P.index]=F;const k=(v,I)=>{const E=h[v],L=C[v];if(!(L.peakValue>1)||L.peakIndex<0)return null;const z=Math.max(0,Math.floor(I-M)),U=Math.min(E.length,Math.ceil(I+M+1));let R=-1/0,O=-1;for(let X=z;X<U;X++)E[X]>R&&(R=E[X],O=X);if(O<0||!(R>1))return null;const B=Math.max(1e-6,L.peakValue),Y=Math.abs(O-I)<=M,Z=R>=B*.25;return!Y||!Z?null:u(E,O,m,v,y,x,b)};let _=F?b?(F.y-m)/x:(F.x-m)/x:P.peakIndex;for(let v=P.index+1;v<f;v++){const I=k(v,_);I&&(T[v]=I,_=b?(I.y-m)/x:(I.x-m)/x)}_=F?b?(F.y-m)/x:(F.x-m)/x:P.peakIndex;for(let v=P.index-1;v>=0;v--){const I=k(v,_);I&&(T[v]=I,_=b?(I.y-m)/x:(I.x-m)/x)}return T.filter(v=>!!v)};if(s){const h=n.map(f=>f.map((p,m)=>m===0?0:Math.abs(p-f[m-1])));return d(h,a,l,t,r,i,!1)}const c=[];for(let h=0;h<l;h++){const f=new Array(a).fill(0);for(let p=1;p<a;p++)f[p]=Math.abs(n[p][h]-n[p-1][h]);c.push(f)}return d(c,l,a,e,i,r,!0)}function ve(n){if(n.length<2)return null;let t=0,e=0,i=0;for(const c of n)t+=c.weight,e+=c.x*c.weight,i+=c.y*c.weight;if(t<=0)return null;e/=t,i/=t;let r=0,s=0,o=0;for(const c of n){const g=c.x-e,h=c.y-i;r+=c.weight*g*g,s+=c.weight*h*h,o+=c.weight*g*h}r/=t,s/=t,o/=t;const a=.5*Math.atan2(2*o,r-s);let l=Math.cos(a),u=Math.sin(a);const d=Math.hypot(l,u);return!Number.isFinite(d)||d<=1e-9?null:(l/=d,u/=d,(l<0||Math.abs(l)<=1e-9&&u<0)&&(l=-l,u=-u),{pointX:e,pointY:i,dirX:l,dirY:u,orientation:Math.abs(l)>=Math.abs(u)?1:2})}function so(n,t){if(n.length!==4||t.length!==4||n.some(i=>i.length!==4))return null;const e=n.map((i,r)=>[...i,t[r]]);for(let i=0;i<4;i++){let r=i,s=Math.abs(e[i][i]);for(let a=i+1;a<4;a++){const l=Math.abs(e[a][i]);l>s&&(s=l,r=a)}if(!(s>1e-12))return null;if(r!==i){const a=e[i];e[i]=e[r],e[r]=a}const o=e[i][i];for(let a=i;a<=4;a++)e[i][a]/=o;for(let a=0;a<4;a++){if(a===i)continue;const l=e[a][i];if(!(Math.abs(l)<=1e-12))for(let u=i;u<=4;u++)e[a][u]-=l*e[i][u]}}return[e[0][4],e[1][4],e[2][4],e[3][4]]}function ao(n){if(n.length<4)return 0;const t=[...n].sort((g,h)=>g.x-h.x),e=t[0].x,r=t[t.length-1].x-e;if(!(r>1e-6))return 0;const s=16,o=[];for(let g=0;g<s;g++){const h=Math.max(0,Math.floor((g-1.5)*t.length/s)),f=Math.min(t.length-1,Math.floor((g+2.5)*t.length/s));if(f<h)continue;let p=0,m=0,y=0;for(let x=h;x<=f;x++)p+=t[x].x,m+=t[x].y,y++;y>0&&o.push({x:p/y,y:m/y})}if(o.length<4)return 0;const a=[.05952381,0,-.03571429,-.04761905,-.03571429,0,.05952381],l=new Array(o.length).fill(0),u=3;for(let g=0;g<o.length;g++){let h=0;for(let f=-u;f<=u;f++){const p=g+f;p<0||p>=o.length||(h+=a[f+u]*o[p].y)}l[g]=h}let d=0,c=1/0;for(let g=0;g<l.length-1;g++){const h=l[g],f=l[g+1];if(h===0){const b=Math.abs(o[g].x);b<c&&(c=b,d=o[g].x);continue}if(h*f>=0)continue;const p=f-h,m=Math.abs(p)>1e-12?-h/p:.5,y=o[g].x+(o[g+1].x-o[g].x)*m,x=Math.abs(y);x<c&&(c=x,d=y)}return!Number.isFinite(d)||d<e+.3*r||d>e+.7*r?0:d}function oo(n){if(n.length<8)return null;const t=[...n].filter(f=>Number.isFinite(f.x)&&Number.isFinite(f.y)&&Number.isFinite(f.weight)).sort((f,p)=>f.x-p.x);if(t.length<8)return null;const i=[ao(t),0,.5*(t[Math.floor((t.length-1)*.5)].x+t[Math.ceil((t.length-1)*.5)].x)];let r=null;for(const f of i){if(!Number.isFinite(f))continue;let p=0,m=0;for(const y of t)y.x<=f?p++:m++;if(p>=4&&m>=4){r=f;break}}if(r===null)return null;const s=Array.from({length:4},()=>new Array(4).fill(0)),o=new Array(4).fill(0);for(const f of t){const p=f.x,m=f.y,y=Math.max(1e-6,f.weight),x=p<=r?[p*p,p,1,0]:[2*r*p-r*r,p,1,(p-r)*(p-r)];for(let b=0;b<4;b++){o[b]+=y*x[b]*m;for(let M=0;M<4;M++)s[b][M]+=y*x[b]*x[M]}}const a=so(s,o);if(!a)return null;const[l,u,d,c]=a,g=u+2*(l-c)*r,h=d+(c-l)*r*r;return[l,u,d,c,g,h].every(f=>Number.isFinite(f))?{splitX:r,left:[l,u,d],right:[c,g,h]}:null}function lo(n,t,e){const[i,r,s]=e;if(Math.abs(i)<=1e-12){const b=1+r*r;return b>1e-12?[(n-r*(s-t))/b]:[n]}const o=2*i*i,a=3*i*r,l=1+2*i*s-2*i*t+r*r,u=r*s-t*r-n;if(Math.abs(o)<=1e-12)return[n];const d=a/o,c=l/o,g=u/o,h=(d*d-3*c)/9,f=(2*d*d*d-9*d*c+27*g)/54,p=f*f-h*h*h;if(p<0&&h>0){const b=Math.acos(Math.max(-1,Math.min(1,f/Math.sqrt(h*h*h)))),M=-2*Math.sqrt(h);return[M*Math.cos(b/3)-d/3,M*Math.cos((b+2*Math.PI)/3)-d/3,M*Math.cos((b-2*Math.PI)/3)-d/3]}const m=Math.sqrt(Math.max(0,p)),y=-Math.sign(f||1)*Math.cbrt(Math.abs(f)+m),x=Math.abs(y)<=1e-12?0:h/y;return[y+x-d/3]}function xs(n,t){if(n.length<8)return null;const e=-t.dirY,i=t.dirX,r=n.map(o=>({x:(o.x-t.pointX)*t.dirX+(o.y-t.pointY)*t.dirY,y:(o.x-t.pointX)*e+(o.y-t.pointY)*i,weight:o.weight})),s=oo(r);return s?{...t,normalX:e,normalY:i,splitX:s.splitX,left:s.left,right:s.right}:null}function co(n,t){const e=n.x-t.pointX,i=n.y-t.pointY,r=e*t.dirX+i*t.dirY,s=e*t.normalX+i*t.normalY,o=r<t.splitX?t.left:t.right,a=lo(r,s,o);let l=s,u=Number.POSITIVE_INFINITY;for(const d of a){if(!Number.isFinite(d))continue;const c=o[0]*d*d+o[1]*d+o[2],g=r-d,h=s-c,f=Math.hypot(g,h);Number.isFinite(f)&&f<u&&(u=f,l=(h>=0?1:-1)*f)}return Number.isFinite(u)?l:s}function bs(n,t,e,i,r,s,o,a=Et){if(!t||t.length<8||n.length===0)return null;const l=t.filter(v=>Number.isFinite(v.x)&&Number.isFinite(v.y)).map(v=>({x:v.x,y:v.y,weight:1}));if(l.length<8)return null;const u=ve(l);if(!u)return null;const d=e.p2.x-e.p1.x,c=e.p2.y-e.p1.y,g=Math.hypot(d,c);if(!Number.isFinite(g)||g<=1e-6)return null;let h=u.dirX,f=u.dirY;h*d+f*c<0&&(h=-h,f=-f);const p={...u,dirX:h,dirY:f},m=xs(l,p),y=d/g,x=c/g,b=-x,M=y,C=(e.p1.x+e.p2.x)*.5,S=(e.p1.y+e.p2.y)*.5,w=-p.dirY,P=p.dirX,T=Math.abs(y)>=Math.abs(x)?1:2,F=[],k=[];for(const v of n){const I=v.x-C,E=v.y-S,L=I*y+E*x;if(Math.abs(L)>i)continue;const z=I*b+E*M;if(Math.abs(z)>r)continue;const U=m?co(v,m):(v.x-p.pointX)*w+(v.y-p.pointY)*P;Number.isFinite(U)&&(F.push(U),k.push(v.value))}if(F.length<8)return null;const _=o?s!=null&&s.forceLegacyModel?kn(F,k,T,a,r*2):yn(F,k,T,a):ci(F,k,Math.max(2,r*2),s==null?void 0:s.manualBinSize,T,s==null?void 0:s.preferAutoPerEdgeBin);return _?(_.quadraticProjectionUsed=!!m,_):null}function Rr(n){if(n.length<2)return null;const t=n.filter(e=>Number.isFinite(e.x)&&Number.isFinite(e.y)).map(e=>({x:e.x,y:e.y,weight:1}));return t.length<2?null:Sn(t,ve(t))}function Sn(n,t){if(!t||n.length<2)return null;let e=1/0,i=-1/0;for(const s of n){const o=(s.x-t.pointX)*t.dirX+(s.y-t.pointY)*t.dirY;e=Math.min(e,o),i=Math.max(i,o)}if(!Number.isFinite(e)||!Number.isFinite(i))return null;const r=Math.max(.5,(i-e)*.03);return{p1:{x:t.pointX+t.dirX*(e-r),y:t.pointY+t.dirY*(e-r)},p2:{x:t.pointX+t.dirX*(i+r),y:t.pointY+t.dirY*(i+r)}}}function uo(n,t,e,i,r,s,o,a){return[Sn(n,t),Sn(e,i),Sn(r,s),Sn(o,a)]}function ho(n,t,e){if(!n||n.length<8)return;const i=n.filter(x=>Number.isFinite(x.x)&&Number.isFinite(x.y)).map(x=>({x:x.x,y:x.y,weight:1}));if(i.length<8)return;const r=ve(i);if(!r)return;const s=t.p2.x-t.p1.x,o=t.p2.y-t.p1.y,a=Math.hypot(s,o);if(!Number.isFinite(a)||a<=1e-6)return;let l=r.dirX,u=r.dirY;l*s+u*o<0&&(l=-l,u=-u);const d=xs(i,{...r,dirX:l,dirY:u});if(!d)return;const c=n.map(x=>(x.x-d.pointX)*d.dirX+(x.y-d.pointY)*d.dirY).filter(x=>Number.isFinite(x)),g=(t.p1.x-d.pointX)*d.dirX+(t.p1.y-d.pointY)*d.dirY,h=(t.p2.x-d.pointX)*d.dirX+(t.p2.y-d.pointY)*d.dirY;if(Number.isFinite(g)&&c.push(g),Number.isFinite(h)&&c.push(h),c.length<2)return;const f=Math.min(...c),p=Math.max(...c);if(!Number.isFinite(f)||!Number.isFinite(p)||p-f<=1e-6)return;const m=Math.max(21,e),y=[];for(let x=0;x<m;x++){const b=m===1?.5:x/(m-1),M=f+(p-f)*b,C=M<d.splitX?d.left:d.right,S=C[0]*M*M+C[1]*M+C[2];y.push({x:d.pointX+M*d.dirX+S*d.normalX,y:d.pointY+M*d.dirY+S*d.normalY})}return y}function ci(n,t,e,i,r,s=!1,o=!1,a=!1){if(n.length===0||t.length!==n.length)return null;const l=fs(n,t),u=l.dists,d=l.vals;if(u.length===0)return null;const c=()=>{const w=e/2;let P=0;for(const F of n)Math.abs(F)<=w&&P++;if(P<=0)return .125;const T=40*w/P;return Math.max(.01,Math.min(.125,T))},g=(w,P,T,F,k)=>{if(!(F>0)||!(k>0)||!(T>P))return!1;const _=Math.floor((T-P)/F);if(_<2)return!1;const v=Math.max(P,-k),I=Math.min(T,k);if(!(I>v))return!1;const E=Math.max(0,Math.floor((v-P)/F)),L=Math.min(_,Math.ceil((I-P)/F));if(L<=E)return!1;const z=new Array(L-E).fill(0),U=P+E*F,R=P+L*F;for(let O=0;O<w.length;O++){const B=w[O];if(B<U)continue;if(B>=R)break;const Y=Math.floor((B-P)/F);Y>=E&&Y<L&&z[Y-E]++}return z.every(O=>O>0)},h=()=>{const w=u[0],P=u[u.length-1],T=Math.max(0,e*.25),F=.125,k=.5,_=.001,v=Math.round((k-F)/_);for(let I=0;I<=v;I++){const E=Number((F+I*_).toFixed(3));if(g(u,w,P,E,T))return E}return k};let f=.125;i&&i>0?f=Math.max(.01,Math.min(.5,i)):s?f=h():f=c();const p=u[0],m=u[u.length-1],y=Math.floor((m-p)/f);if(y<2)return null;const x=()=>{const w=new Array(y).fill(0),P=new Array(y).fill(0);for(let F=0;F<u.length;F++){const k=(u[F]-p)/f;if(Number.isFinite(k))if(o){const _=Math.floor(k),v=k-_,I=1-v,E=v;_>=0&&_<y&&(w[_]+=d[F]*I,P[_]+=I);const L=_+1;L>=0&&L<y&&(w[L]+=d[F]*E,P[L]+=E)}else{const _=Math.floor(k);_>=0&&_<y&&(w[_]+=d[F],P[_]++)}}let T=d[0];for(let F=0;F<y;F++)P[F]>0?(w[F]/=P[F],T=w[F]):w[F]=T;return w},b=a?null:Wa(u,d,f),M=(b==null?void 0:b.esf)??x(),C=(b==null?void 0:b.lsfFull)??Ga(M,f),S=Math.max(0,Math.min(y-1,-p/f-.5));return{esf:M,lsfFull:C,binSize:f,orientation:r,zeroIndex:S,shortSidePx:e,fallbackUsed:a||!b}}function yn(n,t,e,i=Et){if(n.length===0||t.length!==n.length)return null;const r=fs(n,t),s=r.dists.map((a,l)=>({dist:a,value:r.vals[l]})).filter(a=>Math.abs(a.dist)<i);if(s.length<8)return null;const o=us(s.map(a=>a.dist),s.map(a=>a.value),i);return o?{esf:o.esf,lsfFull:o.lsfFull,binSize:li,orientation:e,zeroIndex:se/2,shortSidePx:i*2,fallbackUsed:!1,mtfmapperLike:!0,mtfmapperOrderedDists:s.map(a=>a.dist),mtfmapperOrderedVals:s.map(a=>a.value),mtfmapperEffectiveMaxDot:i}:null}function kn(n,t,e,i=Et,r=i*2){if(n.length===0||t.length!==n.length)return null;const s=[],o=[];for(let a=0;a<n.length;a++){const l=n[a],u=t[a];!Number.isFinite(l)||!Number.isFinite(u)||Math.abs(l)>=i||(s.push(l),o.push(u))}return s.length<8?null:ci(s,o,Math.max(2,r),void 0,e,!0,!0,!0)}function Se(n,t,e,i,r,s,o,a,l=Et){if(s<=0||o<=0)return null;const d=!(!!(a!=null&&a.isThreePlane)&&n.length>=t*e*3)&&((a==null?void 0:a.greenOnly)??!1),c=(a==null?void 0:a.bayerPattern)||"RGGB",g=r.p2.x-r.p1.x,h=r.p2.y-r.p1.y,f=Math.hypot(g,h);if(!Number.isFinite(f)||f<=1e-6)return null;const p=g/f,m=h/f,y=-m,x=p,b=(r.p1.x+r.p2.x)*.5,M=(r.p1.y+r.p2.y)*.5,C=Math.abs(p)>=Math.abs(m)?1:2,S=d?ms(n,t,e,i,0,0,c,a==null?void 0:a.greenPhase,a==null?void 0:a.blackLevel):gs(n,t,e,i,{...a,globalX:0,globalY:0});if(S.length===0)return null;if(!(a!=null&&a.disableQuadraticProjection)){const T=bs(S,a==null?void 0:a.quadraticFitPoints,r,s,o,a,!0,l);if(T)return T}const w=[],P=[];for(const T of S){const F=T.x-b,k=T.y-M,_=F*p+k*m;if(Math.abs(_)>s)continue;const v=F*y+k*x;Math.abs(v)>o||(w.push(v),P.push(T.value))}return w.length<8?null:a!=null&&a.forceLegacyModel?kn(w,P,C,l):yn(w,P,C,l)}function fo(n,t,e=0){const i=[...n.lsfFull];if(i.length<3)return!1;const r=Math.max(n.binSize,1e-6),s=Number.isFinite(n.zeroIndex)?n.zeroIndex:i.length/2,o=Math.max(1,Math.round((n.shortSidePx??0)*.5/r));let{peakPos:a,peakIdx:l,peakVal:u}=ti(i,s,o);const d=u*.2;let c=0,g=i.length-1;for(let p=l;p>=0;p--)if(i[p]<d){c=p;break}for(let p=l;p<i.length;p++)if(i[p]<d){g=p;break}const h=g-c;if(t&&h>0){const p=h*4,m=[],y=[];if(e>0){const x=Math.max(0,l-p-e),b=Math.max(0,l-p);for(let S=x;S<b;S++)m.push(S),y.push(i[S]);const M=Math.min(i.length,l+p),C=Math.min(i.length,l+p+e);for(let S=M;S<C;S++)m.push(S),y.push(i[S])}else{for(let x=0;x<Math.max(0,l-p);x++)m.push(x),y.push(i[x]);for(let x=Math.min(i.length,l+p);x<i.length;x++)m.push(x),y.push(i[x])}if(m.length>2){const{slope:x,intercept:b}=$i(m,y);for(let M=0;M<i.length;M++)i[M]=i[M]-(x*M+b);({peakPos:a}=ti(i,s,o))}}return Math.abs(a-s)*r<=Math.max(1e-6,(n.shortSidePx??0)/6)}function po(n){const t=n.length;if(t<3)return!1;let e=0,i=-1/0;for(let o=0;o<t;o++){const a=Math.abs(n[o]);a>i&&(i=a,e=o)}const r=t/3,s=2*t/3;return e>=r&&e<=s}function _t(n,t,e){const i=Math.max(0,Math.floor(n.x)),r=Math.max(0,Math.floor(n.y)),s=Math.min(t,Math.ceil(n.x+n.w)),o=Math.min(e,Math.ceil(n.y+n.h)),a=s-i,l=o-r;return a<2||l<2?null:{x:i,y:r,w:a,h:l}}function Zi(n,t,e,i){const r=[],s=n.x,o=n.y,a=n.x+n.w,l=n.y+n.h,u=n.x+n.w*.5,d=n.y+n.h*.5,c=[{x:s,y:o},{x:a,y:o},{x:a,y:l},{x:s,y:l},{x:u,y:o},{x:a,y:d},{x:u,y:l},{x:s,y:d},{x:u,y:d}];for(const g of c){const h=Ce(g,t);Number.isFinite(h.x)&&Number.isFinite(h.y)&&r.push(h)}return r.length===0?null:_t(Xt(r,2),e,i)}function Xt(n,t=0){let e=1/0,i=1/0,r=-1/0,s=-1/0;for(const o of n)e=Math.min(e,o.x),i=Math.min(i,o.y),r=Math.max(r,o.x),s=Math.max(s,o.y);return{x:e-t,y:i-t,w:r-e+t*2,h:s-i+t*2}}function Er(n,t){let e=Math.atan2(t,n)*180/Math.PI;return e<0&&(e+=180),e}function xe(n,t){const e=n.p2.x-n.p1.x,i=n.p2.y-n.p1.y,r=Math.hypot(e,i);if(!Number.isFinite(r)||r<=1e-6)return null;const s=-i/r,o=e/r;return[{x:n.p1.x+s*t,y:n.p1.y+o*t},{x:n.p2.x+s*t,y:n.p2.y+o*t},{x:n.p2.x-s*t,y:n.p2.y-o*t},{x:n.p1.x-s*t,y:n.p1.y-o*t}]}function mo(n,t,e,i,r,s,o,a){if(s<=0||o<=0)return null;const u=!(!!(a!=null&&a.isThreePlane)&&n.length>=t*e*3)&&((a==null?void 0:a.greenOnly)??!1),d=(a==null?void 0:a.bayerPattern)||"RGGB",c=r.p2.x-r.p1.x,g=r.p2.y-r.p1.y,h=Math.hypot(c,g);if(!Number.isFinite(h)||h<=1e-6)return null;const f=c/h,p=g/h,m=-p,y=f,x=(r.p1.x+r.p2.x)*.5,b=(r.p1.y+r.p2.y)*.5,M=Math.abs(f)>=Math.abs(p)?1:2,C=u?ms(n,t,e,i,0,0,d,a==null?void 0:a.greenPhase,a==null?void 0:a.blackLevel):gs(n,t,e,i,{...a,globalX:0,globalY:0});if(C.length===0)return null;if(!(a!=null&&a.disableQuadraticProjection)){const P=bs(C,a==null?void 0:a.quadraticFitPoints,r,s,o,a,!1);if(P)return P}const S=[],w=[];for(const P of C){const T=P.x-x,F=P.y-b,k=T*f+F*p;if(Math.abs(k)>s)continue;const _=T*m+F*y;Math.abs(_)>o||(S.push(_),w.push(P.value))}return S.length<8?null:ci(S,w,o*2,a==null?void 0:a.manualBinSize,M,a==null?void 0:a.preferAutoPerEdgeBin)}function go(n,t,e){const i=[...n],r=new Array(n.length).fill(0),s=[0,0,0];let o=-1,a=1,l=-1;for(let c=1;c<n.length-1;c++){let g=0;if(n[c]>1e-4){g=Math.atan2(e[c]*o,t[c]*o);let h=0;for(let f=-5;f<=5;f++)Math.abs(g+f*2*Math.PI-s[1])<Math.abs(g+h*2*Math.PI-s[1])&&(h=f);g+=h*2*Math.PI}c>3&&Math.abs(g-s[0])>Math.PI/2&&l<c-1&&n[c]<.5&&(a*=-1,l=c),i[c]*=a,o*=-1,s[0]=s[1],s[1]=g,s[2]=g}const u=[-.086,.343,.486,.343,-.086];for(let c=0;c<n.length-3;c++){let g=0;for(let h=-2;h<=2;h++)g+=i[Math.abs(c+h)]*u[h+2];r[c]=g}for(let c=0;c<n.length-3;c++)i[c]=r[c];const d=7;for(let c=0;c<3;c++){r.fill(0);for(let h=0;h<n.length-d;h++)if(h<d)r[h]=i[h];else{const f=Math.min(5,Math.floor((h-5)/3)),p=Ya[f];let m=0;for(let y=-d;y<=d;y++)m+=i[h+y]*p[y+d];r[h]=m}for(let h=n.length-d-2;h<n.length;h++)r[h]=i[h];const g=Math.abs(r[0])>1e-9?r[0]:1;for(let h=0;h<n.length;h++)i[h]=r[h]/g}for(let c=0;c<n.length;c++)i[c]=Math.abs(i[c]);return i}function yo(n,t){const e=[[0,0,0],[0,0,0],[0,0,0]],i=[0,0,0];for(let o=0;o<n.length;o++){const a=n[o],l=-t+o,u=[1,a,a*a];for(let d=0;d<3;d++){i[d]+=u[d]*l;for(let c=0;c<3;c++)e[d][c]+=u[d]*u[c]}}const r=Es(e);if(!r)return null;const s=Ls(r,i);return[s[0],s[1],s[2]]}function xo(n,t){let e=0,i=1,r=0,s=!1,o=0;const a=Math.min(n.length,se/16*2);for(let l=0;l<a&&!s;l++){const u=n[l];if(i>.5&&u<=.5){const d=-(u-i)*se;Math.abs(d)>1e-9&&(r=-((.5-i-d*e)/d),o=l,s=!0)}i=u,e=l/se}if(!s)return null;if(o>=5&&o<a-10){const l=Math.min(Math.max(2,o-9),9),u=yo(n.slice(o-l,o+l+1),l);if(u){const c=(u[0]+.5*u[1]+.25*u[2]+o)/se;if(o>9)r=c;else{const h=(o-5)/se/8;r=(1-h)*r+h*c}}}return r*Lt*t}function bo(n,t){if(n.length===0)return null;const e=se,i=se/16*4,r=new Cn(e),s=1,o=ja(),a=new Float32Array(501);for(let _=0;_<=500;_++)a[_]=_/500*s*2;const l=new Array(i).fill(0).map((_,v)=>v/e*s*Lt),u=new Float32Array(i).fill(0),d=new Float32Array(i).fill(0);let c=0,g=[],h=[],f=[],p=[],m=[],y=[],x=[],b=null,M=0;for(const _ of n){const v=_.mtfmapperOrderedDists&&_.mtfmapperOrderedVals&&_.mtfmapperOrderedDists.length===_.mtfmapperOrderedVals.length?us(_.mtfmapperOrderedDists,_.mtfmapperOrderedVals,_.mtfmapperEffectiveMaxDot??Et):null,I=(v==null?void 0:v.lsfFull)??_.lsfFull,E=(v==null?void 0:v.esf)??_.esf;if(I.length<e)continue;const L=new Float32Array(e);for(let R=0;R<e;R++)L[R]=I[R]??0;r.transform(L);const z=Math.max(1e-9,Math.abs(r._real[0])),U=new Array(i).fill(0);for(let R=1;R<i;R++)U[R]=Math.atan2(r._imag[R],r._real[R]);for(let R=0;R<i;R++)u[R]+=r._real[R]/z,d[R]+=r._imag[R]/z;if(c++,M+=_.shortSidePx*.5,g.length===0){g=[...I],h=[...E];const R=new Array(i).fill(0);R[0]=1;for(let X=1;X<i;X++)R[X]=Math.hypot(r._real[X]/z,r._imag[X]/z);const O=l.map(X=>X),B=(Number.isFinite(_.zeroIndex)?_.zeroIndex:0)*(_.binSize??li),Y=hs(U,O,R,Number.POSITIVE_INFINITY,B),Z=ds(Y,l,a);p=Z.ptfRaw,m=Z.ptfUnwrapped,y=Z.ptfLinear,x=Z.ptfResidual,f=Z.ptfResidual,b=Y.fit}}if(c===0)return null;const C=new Float32Array(i),S=new Float32Array(i),w=new Array(i).fill(0);w[0]=1;for(let _=0;_<i;_++)C[_]=u[_]/c,S[_]=d[_]/c,_>0&&(w[_]=Math.hypot(C[_],S[_]));const P=go(w,C,S),T=new Array(i).fill(0);for(let _=0;_<i;_++)T[_]=P[_]/o[_];const F=Array.from(a,_=>Je(_,l,T)),k=xo(T,s);return{esf:h,lsf:[],lsfCropped:g,mtf:F,ptf:f,ptfRaw:p,ptfUnwrapped:m,ptfLinear:y,ptfResidual:x,ptfPhaseFit:b,freqs:Array.from(a),mtf50:k,calcRadius:M/c}}function Mo(n,t,e,i=!1,r=0,s=!1){if(n.length===0)return null;if(n.every(w=>w.mtfmapperLike))return bo(n);const o=4096,a=new Cn(o),l=1,u=new Float32Array(501);for(let w=0;w<=500;w++)u[w]=w/500*l*2;const d=new Float32Array(501).fill(0);let c=0,g=[],h=[],f=0,p=[],m=[],y=[],x=[],b=[],M=null;for(const w of n){let P=[...w.lsfFull];const T=w.binSize,F=Number.isFinite(w.zeroIndex)?w.zeroIndex:P.length/2,k=Math.max(1,Math.round((w.shortSidePx??0)*.5/Math.max(T,1e-6)));let{peakPos:_,peakIdx:v,peakVal:I}=ti(P,F,k);const E=I*.2;let L=0,z=P.length-1;for(let A=v;A>=0;A--)if(P[A]<E){L=A;break}for(let A=v;A<P.length;A++)if(P[A]<E){z=A;break}const U=z-L;let R=!1;if(i&&U>0){const A=U*4,Q=[],G=[];if(r>0){const V=Math.max(0,v-A-r),rt=Math.max(0,v-A);for(let dt=V;dt<rt;dt++)Q.push(dt),G.push(P[dt]);const st=Math.min(P.length,v+A),it=Math.min(P.length,v+A+r);for(let dt=st;dt<it;dt++)Q.push(dt),G.push(P[dt])}else{for(let V=0;V<Math.max(0,v-A);V++)Q.push(V),G.push(P[V]);for(let V=Math.min(P.length,v+A);V<P.length;V++)Q.push(V),G.push(P[V])}if(Q.length>2){const{slope:V,intercept:rt}=$i(Q,G);for(let st=0;st<P.length;st++)P[st]=P[st]-(V*st+rt);({peakPos:_,peakIdx:v,peakVal:I}=ti(P,F,k)),R=!0}}let O=0,B=0;if(t>0)B=t,O=Math.round(t/T);else{const A=I*.2;let Q=0,G=P.length-1;for(let it=v;it>=0;it--)if(P[it]<A){Q=it;break}for(let it=v;it<P.length;it++)if(P[it]<A){G=it;break}const rt=(G-Q)*T;let st=Math.max(2,rt*8);B=st,O=Math.round(st/T)}f+=B;const Y=Math.max(0,Math.floor(F-O)),Z=Math.min(P.length,Math.ceil(F+O)),X=P.slice(Y,Z);if(X.length<8)continue;const J=new Float32Array(o).fill(0),W=new Array(X.length).fill(0);for(let A=0;A<X.length;A++){let Q=1;s&&(Q=.5*(1-Math.cos(2*Math.PI*A/(X.length-1)))),W[A]=X[A]*Q}const tt=Math.max(0,Math.min(X.length-1,_-Y));for(let A=0;A<o;A++)J[A]=eo(W,o,A+tt);a.transform(J);const lt=[],D=[],$=[];for(let A=0;A<=o/2;A++){const Q=a._real[A],G=a._imag[A],V=Math.sqrt(Q*Q+G*G);lt.push(V),D.push(A/(o*T)*l),$.push(Math.atan2(G,Q))}const K=lt[0];if(K>0){for(let A=0;A<=500;A++){const Q=u[A],V=Va(Q,T);d[A]+=Je(Q,D,lt)/K/V}if(c++,g.length===0){g=no(X,tt,(X.length-1)/2),h=R?_o(P):w.esf;const A=lt.map(rt=>rt/K),Q=D.map(rt=>rt),G=hs($,Q,A,Number.POSITIVE_INFINITY,0),V=ds(G,D,u);m=V.ptfRaw,y=V.ptfUnwrapped,x=V.ptfLinear,b=V.ptfResidual,p=V.ptfResidual,M=G.fit}}}if(c===0)return null;const C=Array.from(d).map(w=>w/c);let S=null;for(let w=0;w<C.length-1;w++)if(C[w]>=.5&&C[w+1]<.5){S=u[w]+(.5-C[w])*(u[w+1]-u[w])/(C[w+1]-C[w]);break}return{esf:h,lsf:[],lsfCropped:g,mtf:C,ptf:p,ptfRaw:m,ptfUnwrapped:y,ptfLinear:x,ptfResidual:b,ptfPhaseFit:M,freqs:Array.from(u),mtf50:S,calcRadius:f/c}}function _o(n){const t=new Array(n.length).fill(0);let e=0;for(let i=0;i<n.length;i++)e+=n[i],t[i]=e;return t}function Je(n,t,e){if(n<=t[0])return e[0];if(n>=t[t.length-1])return e[e.length-1];let i=0;for(;n>t[i+1];)i++;const r=(n-t[i])/(t[i+1]-t[i]);return e[i]+r*(e[i+1]-e[i])}function tr(n){return{...Gi,...n,gradientPercentiles:n!=null&&n.gradientPercentiles&&n.gradientPercentiles.length>0?n.gradientPercentiles:Gi.gradientPercentiles}}function wo(n){return!n||n.length===0?void 0:[Number.isFinite(n[0])?n[0]:0,Number.isFinite(n[1])?n[1]:Number.isFinite(n[0])?n[0]:0,Number.isFinite(n[2])?n[2]:Number.isFinite(n[0])?n[0]:0,Number.isFinite(n[3])?n[3]:Number.isFinite(n[0])?n[0]:0]}function So(n,t){const e=n.width,i=n.height,r=n.data,s=n.bayerPattern||"RGGB",o=wo(n.blackLevels),a=new Float32Array(e*i),l=(M,C)=>M<0||C<0||M>=e||C>=i?null:Math.max(0,r[C*e+M]-gn(o,M,C));let u=1/0,d=-1/0;for(let M=0;M<i;M++){const C=M*e;for(let S=0;S<e;S++){const w=C+S;let P=0;if(vt(S,M,s,t))P=l(S,M)??0;else{const T=[],F=l(S-1,M),k=l(S+1,M),_=l(S,M-1),v=l(S,M+1);if(F!==null&&vt(S-1,M,s,t)&&T.push(F),k!==null&&vt(S+1,M,s,t)&&T.push(k),_!==null&&vt(S,M-1,s,t)&&T.push(_),v!==null&&vt(S,M+1,s,t)&&T.push(v),T.length>0)P=Zn(T);else{const I=[],E=l(S-1,M-1),L=l(S+1,M-1),z=l(S-1,M+1),U=l(S+1,M+1);E!==null&&vt(S-1,M-1,s,t)&&I.push(E),L!==null&&vt(S+1,M-1,s,t)&&I.push(L),z!==null&&vt(S-1,M+1,s,t)&&I.push(z),U!==null&&vt(S+1,M+1,s,t)&&I.push(U),P=Zn(I)}}a[w]=P,P<u&&(u=P),P>d&&(d=P)}}if(!Number.isFinite(u)||!Number.isFinite(d)||d<=u+1e-9)return new Uint8Array(e*i);const c=1024,g=new Uint32Array(c),h=d-u;for(let M=0;M<a.length;M++){const C=Math.max(0,Math.min(1,(a[M]-u)/h)),S=Math.min(c-1,Math.max(0,Math.floor(C*(c-1))));g[S]++}const f=a.length,p=M=>{const C=f*M;let S=0;for(let w=0;w<c;w++)if(S+=g[w],S>=C)return u+w/Math.max(1,c-1)*h;return d},m=p(.01),y=p(.99),x=Math.max(1e-9,y-m),b=new Uint8Array(e*i);for(let M=0;M<a.length;M++){const C=Math.max(0,Math.min(1,(a[M]-m)/x));b[M]=Math.round(C*255)}return b}function vo(n,t,e){const i=new Float32Array(n.length),r=new Float32Array(n.length),s=new Float32Array(n.length);for(let o=1;o<e-1;o++)for(let a=1;a<t-1;a++){const l=o*t+a,u=n[(o-1)*t+(a-1)],d=n[(o-1)*t+a],c=n[(o-1)*t+(a+1)],g=n[o*t+(a-1)],h=n[o*t+(a+1)],f=n[(o+1)*t+(a-1)],p=n[(o+1)*t+a],m=n[(o+1)*t+(a+1)],y=-u-2*g-f+(c+2*h+m),x=-u-2*d-c+(f+2*p+m);i[l]=y,r[l]=x,s[l]=Math.hypot(y,x)}return{gx:i,gy:r,magnitude:s}}function Po(n,t){let e=0,i=0;for(let l=0;l<n.length;l++){const u=n[l];!Number.isFinite(u)||u<=1e-6||(e=Math.max(e,u),i++)}if(i===0||e<=1e-6)return[];const r=1024,s=new Uint32Array(r);for(let l=0;l<n.length;l++){const u=n[l];if(!Number.isFinite(u)||u<=1e-6)continue;const d=Math.max(0,Math.min(1,u/e)),c=Math.min(r-1,Math.floor(d*(r-1)));s[c]++}const o=t&&t.length>0?t:Gi.gradientPercentiles,a=[];for(const l of o){const u=i*l;let d=0;for(let c=0;c<r;c++)if(d+=s[c],d>=u){a.push(c/Math.max(1,r-1)*e);break}}return Array.from(new Set(a.filter(l=>l>0))).sort((l,u)=>u-l)}function Co(n,t){const e=new Uint8Array(n.length);for(let i=0;i<n.length;i++)e[i]=n[i]>=t?1:0;return e}const ko=256*256;function Fo(n,t,e){if(n.length>=ko){const s=Oa.compute(n,t,e);if(s)return{gray:s.blurredGray,gradient:{gx:s.gx,gy:s.gy,magnitude:s.magnitude},backend:"webgl"}}const i=Uo(n,t,e),r=vo(i,t,e);return{gray:i,gradient:r,backend:"cpu"}}function To(n,t,e,i){let r=n;for(let s=0;s<i;s++){const o=new Uint8Array(n.length);for(let a=0;a<e;a++)for(let l=0;l<t;l++){let u=0;for(let d=-1;d<=1&&!u;d++){const c=a+d;if(!(c<0||c>=e))for(let g=-1;g<=1;g++){const h=l+g;if(!(h<0||h>=t)&&r[c*t+h]){u=1;break}}}o[a*t+l]=u}r=o}return r}function Ao(n,t,e){const i=new Int32Array(n.length),r=[];let s=1;for(let o=0;o<n.length;o++){if(!n[o]||i[o]!==0)continue;const a=[o];i[o]=s;let l=0,u=t,d=e,c=0,g=0,h=0,f=!1;for(;l<a.length;){const p=a[l++],m=p%t,y=Math.floor(p/t);h++,u=Math.min(u,m),d=Math.min(d,y),c=Math.max(c,m),g=Math.max(g,y),(m===0||y===0||m===t-1||y===e-1)&&(f=!0);for(let x=-1;x<=1;x++)for(let b=-1;b<=1;b++){if(b===0&&x===0)continue;const M=m+b,C=y+x;if(M<0||C<0||M>=t||C>=e)continue;const S=C*t+M;!n[S]||i[S]!==0||(i[S]=s,a.push(S))}}r.push({label:s,x:u,y:d,w:c-u+1,h:g-d+1,area:h,touchesBorder:f}),s++}return{labels:i,components:r}}function Ms(n,t){const e=Math.hypot(n,t);if(!Number.isFinite(e)||e<=1e-9)return null;let i=n/e,r=t/e;return(i<0||Math.abs(i)<=1e-9&&r<0)&&(i=-i,r=-r),{x:i,y:r}}function Te(n,t){if(n.length===0)return 0;const e=[...n].sort((o,a)=>o.value-a.value),i=e.reduce((o,a)=>o+Math.max(0,a.weight),0);if(i<=0)return e[Math.floor((e.length-1)*t)].value;const r=Math.max(0,Math.min(1,t))*i;let s=0;for(const o of e)if(s+=Math.max(0,o.weight),s>=r)return o.value;return e[e.length-1].value}function Lr(n){const t=n.filter(i=>Number.isFinite(i)).sort((i,r)=>i-r);if(t.length===0)return 0;const e=i=>{if(i.length===1)return i[0];if(i.length===2)return(i[0]+i[1])*.5;const r=Math.ceil(i.length*.5);let s=0,o=1/0;for(let a=0;a+r-1<i.length;a++){const l=i[a+r-1]-i[a];l<o&&(o=l,s=a)}return e(i.slice(s,s+r))};return e(t)}function Io(n,t,e,i,r,s,o){const a=[];for(let l=r.y;l<r.y+r.h;l++)for(let u=r.x;u<r.x+r.w;u++){const d=l*s+u;if(n[d]!==t||!e[d])continue;const c=i.magnitude[d];!Number.isFinite(c)||c<=1e-6||a.push({x:u,y:l,weight:c,gx:i.gx[d],gy:i.gy[d]})}return a}function No(n){let t=0,e=0,i=0,r=0,s=0;for(const l of n){t+=l.weight,e+=l.x*l.weight,i+=l.y*l.weight;const u=Math.hypot(l.gx,l.gy);if(!Number.isFinite(u)||u<=1e-6)continue;const d=-l.gy/u,c=l.gx/u;r+=l.weight*(d*d-c*c),s+=l.weight*(2*d*c)}if(t<=0)return null;e/=t,i/=t;const o=.5*Math.atan2(s,r),a=Ms(Math.cos(o),Math.sin(o));return a?{centerX:e,centerY:i,dirX:a.x,dirY:a.y,orthoX:-a.y,orthoY:a.x}:null}function zn(n,t){let e=0,i=0;const r=-t.dirY,s=t.dirX;for(const o of n){const a=(o.x-t.pointX)*r+(o.y-t.pointY)*s;i+=o.weight*a*a,e+=o.weight}return e<=0?1/0:Math.sqrt(i/e)}function Dr(n,t,e,i,r){const s=Math.max(0,Math.min(t-1,i)),o=Math.max(0,Math.min(e-1,r)),a=Math.floor(s),l=Math.floor(o),u=Math.min(t-1,a+1),d=Math.min(e-1,l+1),c=s-a,g=o-l,h=n[l*t+a],f=n[l*t+u],p=n[d*t+a],m=n[d*t+u],y=h+(f-h)*c,x=p+(m-p)*c;return y+(x-y)*g}function Ro(n,t,e,i,r,s,o,a){const l=Dr(n,t,e,i-s*a,r-o*a);return Dr(n,t,e,i+s*a,r+o*a)-l}function Vn(n,t,e,i,r){const s=Math.max(1e-6,e-t);if(n.length===0||!Number.isFinite(s))return{points:[],coverageRatio:0,centerCoverageRatio:0};const o=Math.max(1.5,Math.min(4,s/18)),a=Math.max(1,Math.ceil(s/o)),l=new Map;for(const f of n){const p=i(f);if(!Number.isFinite(p)||p<t||p>e)continue;const m=Math.max(0,Math.min(a-1,Math.floor((p-t)/o))),y=f.weight/(1+Math.abs(r(f))),x=l.get(m);(!x||y>x.score)&&l.set(m,{point:f,score:y})}const u=Array.from(l.values()).sort((f,p)=>i(f.point)-i(p.point)).map(f=>f.point),d=Math.max(0,Math.floor(a*.3)),c=Math.max(d+1,Math.ceil(a*.7));let g=0;for(let f=d;f<c;f++)l.has(f)&&g++;const h=Math.max(1,c-d);return{points:u,coverageRatio:u.length/a,centerCoverageRatio:g/h}}function jn(n,t){const e=n.dirX*t.dirY-n.dirY*t.dirX;if(!Number.isFinite(e)||Math.abs(e)<=1e-6)return null;const i=t.pointX-n.pointX,r=t.pointY-n.pointY,s=(i*t.dirY-r*t.dirX)/e;return{x:n.pointX+n.dirX*s,y:n.pointY+n.dirY*s}}function Eo(n){if(n.length<3)return 0;let t=0;for(let e=0;e<n.length;e++){const i=n[e],r=n[(e+1)%n.length];t+=i.x*r.y-r.x*i.y}return t*.5}function Lo(n,t,e,i,r,s,o,a,l){const u=Io(i,r.label,s,o,r,t),d=u.map(N=>({x:N.x,y:N.y}));if(u.length<l.minEdgePoints)return{candidate:null,failureStage:"min_edge_points",pointsCount:u.length,strongEdgePoints:d};const c=No(u);if(!c)return{candidate:null,failureStage:"dominant_axes",pointsCount:u.length,strongEdgePoints:d};const g=u.map(N=>{const ot=N.x-c.centerX,ut=N.y-c.centerY;return{...N,u:ot*c.dirX+ut*c.dirY,v:ot*c.orthoX+ut*c.orthoY}}),h={x:c.centerX,y:c.centerY},f=Te(g.map(N=>({value:N.u,weight:N.weight})),l.extentQuantileLow),p=Te(g.map(N=>({value:N.u,weight:N.weight})),l.extentQuantileHigh),m=Te(g.map(N=>({value:N.v,weight:N.weight})),l.extentQuantileLow),y=Te(g.map(N=>({value:N.v,weight:N.weight})),l.extentQuantileHigh),x=Math.max(1e-6,Math.max(Math.abs(f),Math.abs(p))),b=Math.max(1e-6,Math.max(Math.abs(m),Math.abs(y))),M=72,C=360/M,S=Array.from({length:M},()=>[]),w=N=>{let ot=N%360;return ot<0&&(ot+=360),ot},P=(N,ot)=>{const ut=Math.abs(w(N)-w(ot));return Math.min(ut,360-ut)};g.forEach(N=>{const ot=N.u/x,ut=N.v/b,Tt=w(Math.atan2(ut,ot)*180/Math.PI),Vt=Math.hypot(ot,ut),Rt=Math.max(0,Math.min(M-1,Math.floor(Tt/C)));S[Rt].push({point:N,angleDeg:Tt,normRadius:Vt})});const T=S.map(N=>N.length>0?Lr(N.map(ot=>ot.normRadius)):-1/0),F=(N,ot)=>{let ut=-1,Tt=-1/0;for(let yt=0;yt<S.length;yt++){if(S[yt].length===0)continue;const pe=(yt+.5)*C;if(P(pe,N)>45||ot.some(Is=>P(pe,Is)<45))continue;const Fe=T[yt];Fe>Tt&&(Tt=Fe,ut=yt)}let Vt=ut>=0?(ut+.5)*C:N,Rt=ut>=0?S[ut]:g.map(yt=>{const fe=yt.u/x,pe=yt.v/b;return{point:yt,angleDeg:w(Math.atan2(pe,fe)*180/Math.PI),normRadius:Math.hypot(fe,pe)}}).filter(yt=>P(yt.angleDeg,N)<=45&&!ot.some(fe=>P(yt.angleDeg,fe)<45));if(Rt.length===0&&(Rt=g.map(yt=>{const fe=yt.u/x,pe=yt.v/b;return{point:yt,angleDeg:w(Math.atan2(pe,fe)*180/Math.PI),normRadius:Math.hypot(fe,pe)}}).filter(yt=>P(yt.angleDeg,N)<=45),Vt=N),Rt.length===0)return{x:g[0].x,y:g[0].y,u:g[0].u,v:g[0].v,angleDeg:N};const Un=ut>=0?T[ut]:Lr(Rt.map(yt=>yt.normRadius));let _e=0,Mn=0,cr=0,ur=0,hr=0;for(const yt of Rt){const fe=P(yt.angleDeg,N)/45,pe=Math.abs(yt.normRadius-Un),Fe=Math.max(1e-6,yt.point.weight)/(1+fe*2+pe*6);_e+=Fe,Mn+=yt.point.x*Fe,cr+=yt.point.y*Fe,ur+=yt.point.u*Fe,hr+=yt.point.v*Fe}return _e>0?{x:Mn/_e,y:cr/_e,u:ur/_e,v:hr/_e,angleDeg:Vt}:{x:Rt[0].point.x,y:Rt[0].point.y,u:Rt[0].point.u,v:Rt[0].point.v,angleDeg:Rt[0].angleDeg}},k=F(225,[]),_=F(315,[k.angleDeg]),v=F(45,[k.angleDeg,_.angleDeg]),I=F(135,[k.angleDeg,_.angleDeg,v.angleDeg]),E=[{x:k.x,y:k.y},{x:_.x,y:_.y},{x:v.x,y:v.y},{x:I.x,y:I.y}],L=p-f,z=y-m,U=Math.min(L,z),R=Math.max(L,z);if(!Number.isFinite(U)||U<l.minSpanPx||R/Math.max(1,U)>l.maxAspectRatio)return{candidate:null,failureStage:"span_aspect",pointsCount:u.length,minSpan:U,maxSpan:R,axisCentroid:h,axisExtremePoints:E,strongEdgePoints:d};const O=Math.max(l.bandMinPx,Math.min(l.bandMaxPx,U*l.bandScale)),B=Math.max(1,Math.min(3,O*.55)),Y=Math.max(a,0),Z=void 0,X=void 0,J=N=>N.map(ot=>({x:ot.x,y:ot.y,weight:ot.weight})),W=N=>N.map(ot=>({x:ot.x,y:ot.y})),tt=(N,ot,ut)=>N.filter(Tt=>{if(!Number.isFinite(Tt.weight)||Tt.weight<Y)return!1;const Vt=Ro(n,t,e,Tt.x,Tt.y,ot,ut,B);return Number.isFinite(Vt)&&Vt>=l.minPointContrast}),lt=f,D=p,$=m,K=y,A=l.minCoverageRatio,Q=l.minCenterCoverageRatio,G=[],V=[],rt=[],st=[],it=[],dt=(N,ot,ut,Tt,Vt,Rt)=>(ut-N)*(Rt-ot)-(Tt-ot)*(Vt-N),ft=N=>N>1e-6?1:N<-1e-6?-1:0,et=[{u:(k.u+_.u)*.5,v:(k.v+_.v)*.5},{u:(_.u+v.u)*.5,v:(_.v+v.v)*.5},{u:(v.u+I.u)*.5,v:(v.v+I.v)*.5},{u:(I.u+k.u)*.5,v:(I.v+k.v)*.5}],Bt=(N,ot)=>{const ut=ft(dt(k.u,k.v,v.u,v.v,N,ot)),Tt=ft(dt(_.u,_.v,I.u,I.v,N,ot));return`${ut},${Tt}`},Gt=new Map;et.forEach((N,ot)=>{Gt.set(Bt(N.u,N.v),ot)});for(const N of g){if(!Number.isFinite(N.u)||!Number.isFinite(N.v)){it.push(N);continue}let ut=Gt.get(Bt(N.u,N.v))??-1;if(ut<0){let Tt=1/0;for(let Vt=0;Vt<et.length;Vt++){const Rt=et[Vt],Un=(N.u-Rt.u)/x,_e=(N.v-Rt.v)/b,Mn=Un*Un+_e*_e;Mn<Tt&&(Tt=Mn,ut=Vt)}}ut===0?G.push(N):ut===1?V.push(N):ut===2?rt.push(N):ut===3?st.push(N):it.push(N)}const j=[...G,...rt],ct=[...V,...st],at={dir:j.length,ortho:ct.length,unassigned:g.length-j.length-ct.length},pt=G.length>=l.minSidePoints?Te(G.map(N=>({value:N.v,weight:N.weight})),.5):m,mt=rt.length>=l.minSidePoints?Te(rt.map(N=>({value:N.v,weight:N.weight})),.5):y,te=st.length>=l.minSidePoints?Te(st.map(N=>({value:N.u,weight:N.weight})),.5):f,Ut=V.length>=l.minSidePoints?Te(V.map(N=>({value:N.u,weight:N.weight})),.5):p,Mt=[{x:(k.x+_.x)*.5,y:(k.y+_.y)*.5},{x:(_.x+v.x)*.5,y:(_.y+v.y)*.5},{x:(v.x+I.x)*.5,y:(v.y+I.y)*.5},{x:(I.x+k.x)*.5,y:(I.y+k.y)*.5}],le=G.filter(N=>Math.abs(N.v-pt)<=O&&N.u>=lt&&N.u<=D),Kt=rt.filter(N=>Math.abs(N.v-mt)<=O&&N.u>=lt&&N.u<=D),ke=st.filter(N=>Math.abs(N.u-te)<=O&&N.v>=$&&N.v<=K),q=V.filter(N=>Math.abs(N.u-Ut)<=O&&N.v>=$&&N.v<=K),Ft=[le.length,q.length,Kt.length,ke.length],At=[W(le),W(q),W(Kt),W(ke)],xt=tt(le,-c.orthoX,-c.orthoY),wt=tt(Kt,c.orthoX,c.orthoY),ee=tt(ke,-c.dirX,-c.dirY),$t=tt(q,c.dirX,c.dirY),xn=[xt.length,$t.length,wt.length,ee.length],Be=[W(xt),W($t),W(wt),W(ee)],Ue=Vn(xt,lt,D,N=>N.u,N=>N.v-pt),Oe=Vn($t,$,K,N=>N.v,N=>N.u-Ut),Xe=Vn(wt,lt,D,N=>N.u,N=>N.v-mt),Ge=Vn(ee,$,K,N=>N.v,N=>N.u-te),Fn=(N,ot)=>N.slice().sort((ut,Tt)=>ot(ut)-ot(Tt)),Tn=Fn(xt,N=>N.u),An=Fn($t,N=>N.v),In=Fn(wt,N=>N.u),Nn=Fn(ee,N=>N.v),ui=[Tn.length,An.length,In.length,Nn.length],Wt=[Ue.coverageRatio,Oe.coverageRatio,Xe.coverageRatio,Ge.coverageRatio];Ue.centerCoverageRatio,Oe.centerCoverageRatio,Xe.centerCoverageRatio,Ge.centerCoverageRatio;const Yt=[W(Tn),W(An),W(In),W(Nn)],Jt={axisPointCounts:at,sideBandPointCounts:Ft,sideContrastPointCounts:xn,gradientThreshold:a,pointAxisMinDot:Z,pointAxisMargin:X,bandWidth:O,minPointContrast:l.minPointContrast,minCoverageRatio:A,minCenterCoverageRatio:Q,axisCentroid:h,axisExtremePoints:E,axisSideCenters:Mt,strongEdgePoints:d,axisDirPoints:W(j),axisOrthoPoints:W(ct),axisUnassignedPoints:W(it),sideBandPoints:At,sideContrastPoints:Be};if(Tn.length<l.minSidePoints||In.length<l.minSidePoints||Nn.length<l.minSidePoints||An.length<l.minSidePoints)return{candidate:null,failureStage:"min_side_points",pointsCount:u.length,minSpan:U,maxSpan:R,sidePointCounts:ui,sideCoverageRatios:Wt,...Jt,sideFitPoints:Yt};if(Ue.coverageRatio<A||Oe.coverageRatio<A||Xe.coverageRatio<A||Ge.coverageRatio<A||Ue.centerCoverageRatio<Q||Oe.centerCoverageRatio<Q||Xe.centerCoverageRatio<Q||Ge.centerCoverageRatio<Q)return{candidate:null,failureStage:"side_coverage",pointsCount:u.length,minSpan:U,maxSpan:R,sidePointCounts:ui,sideCoverageRatios:Wt,...Jt,sideFitPoints:Yt};const hi=Tn,di=In,fi=Nn,pi=An,ne=ui,Ye=ve(J(hi)),ze=ve(J(di)),Ve=ve(J(fi)),je=ve(J(pi));if(!Ye||!ze||!Ve||!je)return{candidate:null,failureStage:"fit_lines",pointsCount:u.length,minSpan:U,maxSpan:R,sidePointCounts:ne,sideCoverageRatios:Wt,...Jt,sideFitPoints:Yt};const Me=uo(J(hi),Ye,J(pi),je,J(di),ze,J(fi),Ve),Rn=l.minAxisDot,En=(N,ot,ut)=>Math.abs(N.dirX*ot+N.dirY*ut),zt=[En(Ye,c.dirX,c.dirY),En(je,c.orthoX,c.orthoY),En(ze,c.dirX,c.dirY),En(Ve,c.orthoX,c.orthoY)];if(zt[0]<Rn||zt[1]<Rn||zt[2]<Rn||zt[3]<Rn)return{candidate:null,failureStage:"axis_alignment",pointsCount:u.length,minSpan:U,maxSpan:R,sidePointCounts:ne,sideCoverageRatios:Wt,axisDots:zt,...Jt,sideFitPoints:Yt,sideFitLines:Me};const ie=Math.max(l.residualLimitFloor,O*l.residualLimitScale),gt=[zn(J(hi),Ye),zn(J(di),ze),zn(J(fi),Ve),zn(J(pi),je)],mi=[gt[0],gt[3],gt[1],gt[2]],gi=Math.max(...gt),ce=jn(Ye,Ve),ue=jn(Ye,je),he=jn(ze,je),de=jn(ze,Ve);if(!ce||!ue||!he||!de)return{candidate:null,failureStage:"corners",pointsCount:u.length,minSpan:U,maxSpan:R,sidePointCounts:ne,sideCoverageRatios:Wt,axisDots:zt,sideResiduals:[gt[0],gt[3],gt[1],gt[2]],residualLimit:ie,...Jt,sideFitPoints:Yt,sideFitLines:Me};const bn=[ce,ue,he,de],re=Math.abs(Eo(bn));if(!Number.isFinite(re)||re<l.minQuadArea)return{candidate:null,failureStage:"quad_area",pointsCount:u.length,minSpan:U,maxSpan:R,sidePointCounts:ne,sideCoverageRatios:Wt,axisDots:zt,sideResiduals:[gt[0],gt[3],gt[1],gt[2]],residualLimit:ie,...Jt,sideFitPoints:Yt,sideFitLines:Me,quadArea:re};const yi=Math.hypot(ue.x-ce.x,ue.y-ce.y),xi=Math.hypot(he.x-ue.x,he.y-ue.y),bi=Math.hypot(he.x-de.x,he.y-de.y),Mi=Math.hypot(de.x-ce.x,de.y-ce.y),We=[yi,xi,bi,Mi],_i=Math.min(yi,xi,bi,Mi),ks=Math.max(yi,xi,bi,Mi);if(!Number.isFinite(_i)||_i<l.minSideLength||ks/Math.max(1,_i)>l.maxAspectRatio)return{candidate:null,failureStage:"side_length",pointsCount:u.length,minSpan:U,maxSpan:R,sidePointCounts:ne,sideCoverageRatios:Wt,axisDots:zt,sideResiduals:[gt[0],gt[3],gt[1],gt[2]],residualLimit:ie,...Jt,sideFitPoints:Yt,sideFitLines:Me,quadArea:re,sideLengths:We};const He=Ms(ue.x-ce.x+(he.x-de.x),ue.y-ce.y+(he.y-de.y));if(!He)return{candidate:null,failureStage:"corners",pointsCount:u.length,minSpan:U,maxSpan:R,sidePointCounts:ne,sideCoverageRatios:Wt,axisDots:zt,sideResiduals:[gt[0],gt[3],gt[1],gt[2]],residualLimit:ie,...Jt,sideFitPoints:Yt,sideFitLines:Me,quadArea:re,sideLengths:We};const ar={x:-He.y,y:He.x},wi=(ce.x+ue.x+he.x+de.x)*.25,Si=(ce.y+ue.y+he.y+de.y)*.25,Ln=bn.map(N=>{const ot=N.x-wi,ut=N.y-Si;return{u:ot*He.x+ut*He.y,v:ot*ar.x+ut*ar.y}}),Dn=(Math.max(...Ln.map(N=>N.u))-Math.min(...Ln.map(N=>N.u)))*.5,Bn=(Math.max(...Ln.map(N=>N.v))-Math.min(...Ln.map(N=>N.v)))*.5;if(!Number.isFinite(Dn)||!Number.isFinite(Bn)||Math.min(Dn,Bn)<6)return{candidate:null,failureStage:"box_size",pointsCount:u.length,minSpan:U,maxSpan:R,sidePointCounts:ne,sideCoverageRatios:Wt,axisDots:zt,sideResiduals:[gt[0],gt[3],gt[1],gt[2]],residualLimit:ie,sideFitPoints:Yt,quadArea:re,sideLengths:We};const H=Yo(n,t,e,bn,Ye,je,ze,Ve,wi,Si,Dn,Bn,l.innerPurityStdScale,l.outerMeanSpreadLimit);if(!Number.isFinite(gi)||gi>ie)return{candidate:null,failureStage:"residual",pointsCount:u.length,minSpan:U,maxSpan:R,sidePointCounts:ne,sideCoverageRatios:Wt,axisDots:zt,sideResiduals:mi,residualLimit:ie,...Jt,sideFitPoints:Yt,sideFitLines:Me,quadArea:re,sideLengths:We,outerContrast:H.contrast,outerUniformityOk:H.ok,outerMeanSpread:H.meanSpread,outerMeanSpreadLimit:H.meanSpreadLimit,outerAvgStd:H.avgStd,outerAvgStdLimit:H.avgStdLimit,outerSideMeans:H.outerSideMeans,outerSideStds:H.outerSideStds,outerSideStdLimit:H.outerSideStdLimit,outerSideQuads:H.outerSideQuads,innerSideUniformityOk:H.innerSideOk,innerSideStds:H.innerSideStds,innerSideStdLimit:H.innerSideStdLimit,innerSideQuads:H.innerSideQuads};const or=l.filterBlockPurity&&(!H.ok||!H.innerSideOk);if(or||H.contrast<l.minOuterContrast)return{candidate:null,failureStage:or?H.ok?"inner_roi_uniformity":"outer_uniformity":"outer_contrast",pointsCount:u.length,minSpan:U,maxSpan:R,sidePointCounts:ne,sideCoverageRatios:Wt,axisDots:zt,sideResiduals:mi,residualLimit:ie,...Jt,sideFitPoints:Yt,sideFitLines:Me,quadArea:re,sideLengths:We,outerContrast:H.contrast,outerUniformityOk:H.ok,outerMeanSpread:H.meanSpread,outerMeanSpreadLimit:H.meanSpreadLimit,outerAvgStd:H.avgStd,outerAvgStdLimit:H.avgStdLimit,outerSideMeans:H.outerSideMeans,outerSideStds:H.outerSideStds,outerSideStdLimit:H.outerSideStdLimit,outerSideQuads:H.outerSideQuads,innerSideUniformityOk:H.innerSideOk,innerSideStds:H.innerSideStds,innerSideStdLimit:H.innerSideStdLimit,innerSideQuads:H.innerSideQuads};const lr=_t(Xt(bn,1),t,e);if(!lr)return{candidate:null,failureStage:"bbox",pointsCount:u.length,minSpan:U,maxSpan:R,sidePointCounts:ne,sideCoverageRatios:Wt,axisDots:zt,sideResiduals:[gt[0],gt[3],gt[1],gt[2]],residualLimit:ie,...Jt,sideFitPoints:Yt,sideFitLines:Me,quadArea:re,sideLengths:We,outerContrast:H.contrast,outerUniformityOk:H.ok,outerMeanSpread:H.meanSpread,outerMeanSpreadLimit:H.meanSpreadLimit,outerAvgStd:H.avgStd,outerAvgStdLimit:H.avgStdLimit,outerSideMeans:H.outerSideMeans,outerSideStds:H.outerSideStds,outerSideStdLimit:H.outerSideStdLimit,outerSideQuads:H.outerSideQuads,innerSideUniformityOk:H.innerSideOk,innerSideStds:H.innerSideStds,innerSideStdLimit:H.innerSideStdLimit,innerSideQuads:H.innerSideQuads};const Fs=1/(1+gi/Math.max(1,ie)),Ts=l.filterBlockPurity?H.score:1,As=H.contrast*Ts*Fs*Math.sqrt(re);return{candidate:{centerX:wi,centerY:Si,dirX:He.x,dirY:He.y,halfWidth:Dn,halfHeight:Bn,score:As,bbox:lr,corners:bn,sideFitPoints:Yt,outerSideMeans:H.outerSideMeans,outerSideQuads:H.outerSideQuads},failureStage:null,pointsCount:u.length,minSpan:U,maxSpan:R,sidePointCounts:ne,sideCoverageRatios:Wt,axisDots:zt,sideResiduals:mi,residualLimit:ie,...Jt,sideFitPoints:Yt,sideFitLines:Me,quadArea:re,sideLengths:We,outerContrast:H.contrast,outerUniformityOk:H.ok,outerMeanSpread:H.meanSpread,outerMeanSpreadLimit:H.meanSpreadLimit,outerAvgStd:H.avgStd,outerAvgStdLimit:H.avgStdLimit,outerSideMeans:H.outerSideMeans,outerSideStds:H.outerSideStds,outerSideStdLimit:H.outerSideStdLimit,outerSideQuads:H.outerSideQuads,innerSideUniformityOk:H.innerSideOk,innerSideStds:H.innerSideStds,innerSideStdLimit:H.innerSideStdLimit,innerSideQuads:H.innerSideQuads}}function Do(n,t,e,i,r,s,o,a,l){return Lo(n,t,e,i,r,s,o,a,l).candidate}function Bo(n,t,e,i,r,s){const o=tr(r),a=Math.max(i*8,i+128),l=Oo(n,t,e,o.downsampleMaxSide);s==null||s("Detecting candidates: downsampling...",.02),s==null||s("Detecting candidates: edge stage...",.06);const u=Fo(l.gray,l.width,l.height),d=u.gray,c=u.gradient;s==null||s(`Detecting candidates: gradient (${u.backend==="webgl"?"WebGL1":"CPU"})...`,.1);const g=Po(c.magnitude,o.gradientPercentiles),h=l.width*l.height,f=Math.max(o.minComponentAreaPx,Math.round(h*o.minComponentAreaRatio)),p=Math.max(f+1,Math.round(h*o.maxComponentAreaRatio)),m=[],y=Math.max(1,g.reduce((k,_,v)=>k+(v<=1,2),0));let x=0;for(let k=0;k<g.length;k++){const _=g[k],v=Co(c.magnitude,_),I=k<=1?[3,2]:[2,1];for(const E of I){const L=x/y;s==null||s(`Detecting candidates: threshold ${k+1}/${g.length}, dilate ${E}`,.12+.78*L);const z=To(v,l.width,l.height,E),{labels:U,components:R}=Ao(z,l.width,l.height);for(const O of R){if(O.touchesBorder||O.area<f||O.area>p)continue;const B=Do(d,l.width,l.height,U,O,v,c,_,o);if(!B)continue;const Y=1/l.scale,Z=B.corners.map(X=>({x:X.x*Y,y:X.y*Y}));m.push({centerX:B.centerX*Y,centerY:B.centerY*Y,dirX:B.dirX,dirY:B.dirY,halfWidth:B.halfWidth*Y,halfHeight:B.halfHeight*Y,score:B.score,bbox:{x:B.bbox.x*Y,y:B.bbox.y*Y,w:B.bbox.w*Y,h:B.bbox.h*Y},corners:Z,sideFitPoints:B.sideFitPoints?[B.sideFitPoints[0].map(X=>({x:X.x*Y,y:X.y*Y})),B.sideFitPoints[1].map(X=>({x:X.x*Y,y:X.y*Y})),B.sideFitPoints[2].map(X=>({x:X.x*Y,y:X.y*Y})),B.sideFitPoints[3].map(X=>({x:X.x*Y,y:X.y*Y}))]:void 0,outerSideMeans:B.outerSideMeans,outerSideQuads:B.outerSideQuads?[B.outerSideQuads[0].map(X=>({x:X.x*Y,y:X.y*Y})),B.outerSideQuads[1].map(X=>({x:X.x*Y,y:X.y*Y})),B.outerSideQuads[2].map(X=>({x:X.x*Y,y:X.y*Y})),B.outerSideQuads[3].map(X=>({x:X.x*Y,y:X.y*Y}))]:void 0})}m.length>a&&(m.sort((O,B)=>B.score-O.score),m.length=a),x++}}console.log(`[SFR Auto Detect] Candidate pool before dedupe: ${m.length}`),s==null||s(`Detecting candidates: deduplicating (0/${Math.max(1,Math.min(m.length,Math.max(i*4,i+32)))})...`,.94),m.sort((k,_)=>_.score-k.score);const b=Math.max(i*4,i+32),M=m.length>b?m.slice(0,b):m,C=[];if(M.length<=256){console.log(`[SFR Auto Detect] Using simple dedupe for ${M.length} candidates`);for(let k=0;k<M.length;k++){const _=M[k];console.log(`[SFR Auto Detect] Simple dedupe candidate ${k+1}/${M.length}`,_.bbox);const v=M.length<=0?1:k/M.length;if(s==null||s(`Detecting candidates: deduplicating (${k}/${M.length})...`,.94+.05*Math.min(1,v)),!C.some(E=>{const L=Math.hypot(_.centerX-E.centerX,_.centerY-E.centerY),z=Math.max(Math.hypot(_.bbox.w,_.bbox.h),Math.hypot(E.bbox.w,E.bbox.h));return Br(_.bbox,E.bbox)>.28||L<z*.18})&&(C.push(_),C.length>=i))break}return s==null||s("Detecting candidates: deduplicating...",1),C}const S=Math.max(32,Math.round(Math.sqrt(Math.max(1,t*e)/4096))),w=new Map,P=new Set,T=k=>Math.floor(k/S),F=(k,_)=>{if(!Number.isFinite(k.bbox.x)||!Number.isFinite(k.bbox.y)||!Number.isFinite(k.bbox.w)||!Number.isFinite(k.bbox.h)||k.bbox.w<=0||k.bbox.h<=0||k.bbox.w>t*4||k.bbox.h>e*4)return;const v=T(k.bbox.x),I=T(k.bbox.x+k.bbox.w),E=T(k.bbox.y),L=T(k.bbox.y+k.bbox.h);for(let z=E;z<=L;z++)for(let U=v;U<=I;U++){const R=`${U},${z}`,O=w.get(R);O?O.push(_):w.set(R,[_])}};for(let k=0;k<M.length;k++){const _=M[k];if(k===0||k%200===0){const U=M.length<=0?1:k/M.length;s==null||s(`Detecting candidates: deduplicating (${k}/${M.length})...`,.94+.05*Math.min(1,U))}P.clear();const v=T(_.bbox.x),I=T(_.bbox.x+_.bbox.w),E=T(_.bbox.y),L=T(_.bbox.y+_.bbox.h);let z=!1;for(let U=E-1;U<=L+1&&!z;U++)for(let R=v-1;R<=I+1&&!z;R++){const O=w.get(`${R},${U}`);if(O)for(const B of O){if(P.has(B))continue;P.add(B);const Y=C[B];if(!Y)continue;const Z=Math.hypot(_.centerX-Y.centerX,_.centerY-Y.centerY),X=Math.max(Math.hypot(_.bbox.w,_.bbox.h),Math.hypot(Y.bbox.w,Y.bbox.h));if(Br(_.bbox,Y.bbox)>.28||Z<X*.18){z=!0;break}}}if(!z){const U=C.length;if(C.push(_),F(_,U),C.length>=i)break}}return s==null||s("Detecting candidates: deduplicating...",1),C}function Uo(n,t,e){const i=new Uint8Array(n.length);for(let r=0;r<e;r++)for(let s=0;s<t;s++){let o=0,a=0;for(let l=-1;l<=1;l++){const u=r+l;if(!(u<0||u>=e))for(let d=-1;d<=1;d++){const c=s+d;c<0||c>=t||(o+=n[u*t+c],a++)}}i[r*t+s]=Math.round(o/Math.max(1,a))}return i}function Oo(n,t,e,i){const r=Math.max(t,e);if(r<=i)return{gray:n,width:t,height:e,scale:1};const s=i/r,o=Math.max(1,Math.round(t*s)),a=Math.max(1,Math.round(e*s)),l=new Uint8Array(o*a);for(let u=0;u<a;u++){const d=Math.min(e-1,Math.floor(u/s));for(let c=0;c<o;c++){const g=Math.min(t-1,Math.floor(c/s));l[u*o+c]=n[d*t+g]}}return{gray:l,width:o,height:a,scale:s}}function Br(n,t){const e=Math.max(n.x,t.x),i=Math.max(n.y,t.y),r=Math.min(n.x+n.w,t.x+t.w),s=Math.min(n.y+n.h,t.y+t.h),o=Math.max(0,r-e),a=Math.max(0,s-i),l=o*a;if(l<=0)return 0;const u=n.w*n.h+t.w*t.h-l;return u>0?l/u:0}function Ur(n){const t=n.length;if(t===0)return{count:0,mean:0,std:1/0};let e=0;for(const s of n)e+=s;const i=e/t;let r=0;for(const s of n){const o=s-i;r+=o*o}return r/=t,{count:t,mean:i,std:Math.sqrt(Math.max(0,r))}}function Xo(n,t,e,i){return{p1:{x:n.x-t*i,y:n.y-e*i},p2:{x:n.x+t*i,y:n.y+e*i}}}function Or(n,t,e,i,r){return[{x:n.p1.x+t*i,y:n.p1.y+e*i},{x:n.p2.x+t*i,y:n.p2.y+e*i},{x:n.p2.x+t*r,y:n.p2.y+e*r},{x:n.p1.x+t*r,y:n.p1.y+e*r}]}function Go(n,t,e){let i=0;for(let r=0;r<4;r++){const s=e[r],o=e[(r+1)%4],a=(o.x-s.x)*(t-s.y)-(o.y-s.y)*(n-s.x);if(Math.abs(a)<=1e-6)continue;const l=a>0?1:-1;if(i===0)i=l;else if(i!==l)return!1}return!0}function Xr(n,t,e,i){const r=_t(Xt(i,1),t,e);if(!r)return[];const s=[];for(let o=r.y;o<r.y+r.h;o++)for(let a=r.x;a<r.x+r.w;a++)Go(a,o,i)&&s.push(n[o*t+a]);return s}function Yo(n,t,e,i,r,s,o,a,l,u,d,c,g,h){const f=d*2,p=c*2,m=Math.hypot(i[1].x-i[0].x,i[1].y-i[0].y),y=Math.hypot(i[2].x-i[1].x,i[2].y-i[1].y),x=Math.hypot(i[2].x-i[3].x,i[2].y-i[3].y),b=Math.hypot(i[3].x-i[0].x,i[3].y-i[0].y),C=Math.max(...[m,y,x,b]),S=Math.max(2,Math.min(f,p)),w=Math.max(4,C*.25),P=Math.max(2,Math.min(12,S*.22)),T=Math.max(1,Math.min(P,Math.max(1,S*.5-1))),F=1,k=Math.max(8,Math.round(Math.min(w,P*3))),_=[[i[0],i[1],i[1],i[0]],[i[1],i[2],i[2],i[1]],[i[2],i[3],i[3],i[2]],[i[3],i[0],i[0],i[3]]],v=[[i[0],i[1],i[1],i[0]],[i[1],i[2],i[2],i[1]],[i[2],i[3],i[3],i[2]],[i[3],i[0],i[0],i[3]]],I=[],E=[],L=[{corners:[i[0],i[1]],seedLine:r,sideLength:m},{corners:[i[1],i[2]],seedLine:s,sideLength:y},{corners:[i[2],i[3]],seedLine:o,sideLength:x},{corners:[i[3],i[0]],seedLine:a,sideLength:b}];for(let G=0;G<L.length;G++){const V=L[G],rt=Math.max(1,V.sideLength*.5-1),st=Math.max(1,Math.min(rt,w*.5)),it={x:(V.corners[0].x+V.corners[1].x)*.5,y:(V.corners[0].y+V.corners[1].y)*.5},dt=Xo(it,V.seedLine.dirX,V.seedLine.dirY,st),ft=ge(n,t,e,dt.p1,dt.p2,st,Math.max(4,F+Math.max(P,T)+2)),et=(ft==null?void 0:ft.line)||dt,Bt=et.p2.x-et.p1.x,Gt=et.p2.y-et.p1.y,j=Math.hypot(Bt,Gt);if(!Number.isFinite(j)||j<=1e-6)return{ok:!1,score:0,meanSpread:1/0,meanSpreadLimit:1/0,avgStd:1/0,avgStdLimit:1/0,contrast:0,outerMean:0,outerSideMeans:[0,0,0,0],outerSideStds:[1/0,1/0,1/0,1/0],outerSideStdLimit:1/0,outerSideQuads:_,innerSideOk:!1,innerSideStds:[1/0,1/0,1/0,1/0],innerSideStdLimit:1/0,innerSideQuads:v};let ct=-Gt/j,at=Bt/j;const pt={x:(et.p1.x+et.p2.x)*.5,y:(et.p1.y+et.p2.y)*.5};(pt.x-l)*ct+(pt.y-u)*at<0&&(ct=-ct,at=-at);const mt=Or(et,ct,at,F,F+P),te=Or(et,ct,at,-F,-(F+T));_[G]=mt,v[G]=te,I.push(Xr(n,t,e,mt)),E.push(Xr(n,t,e,te))}const z=I.map(Ur);if(z.some(G=>G.count<k||!Number.isFinite(G.std)))return{ok:!1,score:0,meanSpread:1/0,meanSpreadLimit:1/0,avgStd:1/0,avgStdLimit:1/0,contrast:0,outerMean:0,outerSideMeans:[0,0,0,0],outerSideStds:[1/0,1/0,1/0,1/0],outerSideStdLimit:1/0,outerSideQuads:_,innerSideOk:!1,innerSideStds:[1/0,1/0,1/0,1/0],innerSideStdLimit:1/0,innerSideQuads:v};const U=E.map(Ur);if(U.some(G=>G.count<k||!Number.isFinite(G.std)||!Number.isFinite(G.mean)))return{ok:!1,score:0,meanSpread:1/0,meanSpreadLimit:1/0,avgStd:1/0,avgStdLimit:1/0,contrast:0,outerMean:0,outerSideMeans:[0,0,0,0],outerSideStds:[1/0,1/0,1/0,1/0],outerSideStdLimit:1/0,outerSideQuads:_,innerSideOk:!1,innerSideStds:[1/0,1/0,1/0,1/0],innerSideStdLimit:1/0,innerSideQuads:v};const R=z.map(G=>G.mean),O=R.reduce((G,V)=>G+V,0)/R.length,B=U.reduce((G,V)=>G+V.mean,0)/U.length,Y=Math.abs(B-O),Z=Math.max(...R)-Math.min(...R),X=z.reduce((G,V)=>G+V.std,0)/z.length,J=Math.max(0,h),W=Math.max(6,Math.min(20,Y*.45)),tt=R,lt=z.map(G=>G.std),D=Math.max(W,Math.min(30,W*g)),$=U.map(G=>G.std),K=$.every(G=>G<=D),A=Z<=J&&X<=W,Q=1/(1+Z/Math.max(1,J)+X/Math.max(1,W));return{ok:A,score:Q,meanSpread:Z,meanSpreadLimit:J,avgStd:X,avgStdLimit:W,contrast:Y,outerMean:O,outerSideMeans:tt,outerSideStds:lt,outerSideStdLimit:W,outerSideQuads:_,innerSideOk:K,innerSideStds:$,innerSideStdLimit:D,innerSideQuads:v}}function ge(n,t,e,i,r,s,o){const a=r.x-i.x,l=r.y-i.y,u=Math.hypot(a,l);if(!Number.isFinite(u)||u<=1e-6)return null;const d=a/u,c=l/u,g=-c,h=d,f=(i.x+r.x)*.5,p=(i.y+r.y)*.5,m=xe({p1:i,p2:r},o+2);if(!_t(Xt(m||[i,r],2),t,e))return null;const x=Math.max(8,Math.round(s*2)+1),b=Math.max(8,Math.round(o*2)+1),M=x>1?s*2/(x-1):0,C=b>1?o*2/(b-1):0,S=Array.from({length:x},()=>new Array(b).fill(0));for(let _=0;_<x;_++){const v=-s+M*_;for(let I=0;I<b;I++){const E=-o+C*I,L=f+v*d+E*g,z=p+v*c+E*h;S[_][I]=zo(n,t,e,L,z)}}const w=ys(S,-o,-s,C,M,!0);if(w.length<8)return null;const P=w.map(_=>{const v=_.x,I=_.y;return{x:f+I*d+v*g,y:p+I*c+v*h,weight:_.weight}}),T=ve(P);if(!T)return null;let F=T.dirX,k=T.dirY;return F*d+k*c<0&&(F=-F,k=-k),{line:{p1:{x:T.pointX-F*s,y:T.pointY-k*s},p2:{x:T.pointX+F*s,y:T.pointY+k*s}},fitPoints:P.map(_=>({x:_.x,y:_.y}))}}function zo(n,t,e,i,r){if(t<=0||e<=0||n.length!==t*e)return 0;const s=Math.max(0,Math.min(t-1,i)),o=Math.max(0,Math.min(e-1,r)),a=Math.floor(s),l=Math.floor(o),u=Math.min(t-1,a+1),d=Math.min(e-1,l+1),c=s-a,g=o-l,h=n[l*t+a],f=n[l*t+u],p=n[d*t+a],m=n[d*t+u],y=h*(1-c)+f*c,x=p*(1-c)+m*c;return y*(1-g)+x*g}function Vo(n,t,e,i,r){if(i<=0||r<=0||i>=t-1||r>=e-1)return{gx:0,gy:0};const s=r*t+i;return{gx:(n[s+1]-n[s-1])*.5,gy:(n[s+t]-n[s-t])*.5}}function jo(n){if(n.length<20)return null;const t=n.map(w=>Math.max(0,w.weight));let e=0;for(const w of t)e=Math.max(e,w);if(!(e>0))return null;for(let w=0;w<t.length;w++)t[w]/=e;const i=w=>{let P=0,T=0,F=0;for(let W=0;W<n.length;W++){const tt=w[W];tt>0&&(P+=tt,T+=n[W].x*tt,F+=n[W].y*tt)}if(!(P>0))return null;T/=P,F/=P;let k=0,_=0,v=0;for(let W=0;W<n.length;W++){const tt=w[W];if(!(tt>0))continue;const lt=n[W].x-T,D=n[W].y-F;k+=tt*lt*lt,_+=tt*lt*D,v+=tt*D*D}k/=P,_/=P,v/=P;const I=k+v,E=k*v-_*_,L=-I,z=E,U=Math.max(0,L*L-4*z),R=-.5*(L+(L>=0?1:-1)*Math.sqrt(U)),O=Math.abs(R)>1e-12?R:0,B=Math.abs(R)>1e-12?z/R:I,Y=Math.max(O,B);let Z=0,X=1;Math.abs(_)>1e-10?(Z=Y-v,X=_):k>v&&(Z=1,X=0);const J=Math.atan2(-Z,X);return{centroid:{x:T,y:F},angle:J,totalWeight:P}},r=i(t);if(!r)return null;const s=Math.cos(r.angle),o=Math.sin(r.angle),a=new Array(2*16*8).fill(0);for(let w=0;w<n.length;w++){const P=n[w].x-r.centroid.x,T=n[w].y-r.centroid.y,F=P*s+T*o,k=Math.round(F*8+16*8);if(k>=3&&k<a.length-3)for(let _=-3;_<=3;_++)a[k+_]+=t[w]}let l=16*8;for(let w=-5*8+16*8;w<=5*8+16*8;w++)a[w]>a[l]&&(l=w);let u=l-1;for(;u>1&&a[u]>.05*a[l];)u--;let d=l+1;for(;d<a.length-1&&a[d]>.05*a[l];)d++;let c=Math.max(1,u-8);for(;c>1&&a[c]<=a[u];)c--;let g=Math.min(a.length-1,d+8);for(;g<a.length-1&&a[g]<=a[d];)g++;const h=a.slice();for(let w=1;w<h.length;w++)h[w]+=h[w-1];const f=h[h.length-1];if(!(f>0))return null;let p=0;for(let w=1;w<h.length;w++)Math.abs(h[w]-.1*f)<Math.abs(h[p]-.1*f)&&(p=w);let m=h.length-1;for(let w=h.length-2;w>0;w--)Math.abs(h[w]-.9*f)<Math.abs(h[m]-.9*f)&&(m=w);let y=p/8-16,x=m/8-16;const b=x-y;y-=b*.7,x+=b*.7,y=Math.max((c+u)/16-16,y),x=Math.min((g+d)/16-16,x);const M=t.slice();for(let w=0;w<n.length;w++){const P=n[w].x-r.centroid.x,T=n[w].y-r.centroid.y,F=P*s+T*o;M[w]=F>=y&&F<=x?t[w]**4*(1/(10+Math.abs(F))):0}const C=i(M);if(!C)return null;const S=[];for(let w=0;w<n.length;w++)M[w]>0&&S.push({x:n[w].x,y:n[w].y,weight:M[w]});return S.length<8?null:{centroid:C.centroid,angle:C.angle,keptSamples:S}}function Wo(n,t,e,i,r,s=Et){var P;const o=r.x-i.x,a=r.y-i.y,l=Math.hypot(o,a);if(!Number.isFinite(l)||l<=12)return null;const u=o/l,d=a/l,c=5,g=4*s+.5,h=(T,F,k,_,v)=>{const I={x:T.x-F*l*.5,y:T.y-k*l*.5},E={p1:I,p2:{x:I.x+F*l,y:I.y+k*l}},L=xe(E,g+2),z=_t(Xt(L??[E.p1,E.p2],3),t,e),U=[],R=new Map;if(!z)return{reduced:null,scanlines:R};for(let O=z.y;O<z.y+z.h;O++)for(let B=z.x;B<z.x+z.w;B++){const Y=B,Z=O,X=Y-I.x,J=Z-I.y,W=X*F+J*k;if(!(W>c&&W<l-c))continue;const tt=Y-T.x,lt=Z-T.y,D=tt*_+lt*v;if(Math.abs(D)<12){const{gx:$,gy:K}=Vo(n,t,e,B,O),A=$*$+K*K;A>0&&U.push({x:Y,y:Z,weight:A})}if(Math.abs(D)<g){const $=R.get(O);$?(B<$.start&&($.start=B),B>$.end&&($.end=B)):R.set(O,{start:B,end:B})}}return{reduced:jo(U),scanlines:R}};let f={x:(i.x+r.x)*.5,y:(i.y+r.y)*.5},p=u,m=d,y=-m,x=p,b=h(f,p,m,y,x);if(!b.reduced)return null;f=b.reduced.centroid,y=Math.cos(b.reduced.angle),x=Math.sin(b.reduced.angle),p=-x,m=y,p*u+m*d<0&&(p=-p,m=-m,y=-y,x=-x);let M=h(f,p,m,y,x);if(!M.reduced)return null;const C=Math.hypot(M.reduced.centroid.x-f.x,M.reduced.centroid.y-f.y);f=M.reduced.centroid,y=Math.cos(M.reduced.angle),x=Math.sin(M.reduced.angle),p=-x,m=y,p*u+m*d<0&&(p=-p,m=-m,y=-y,x=-x);let S=M;if(C>1){const T=h(f,p,m,y,x);T.reduced&&(S=T,f=T.reduced.centroid,y=Math.cos(T.reduced.angle),x=Math.sin(T.reduced.angle),p=-x,m=y,p*u+m*d<0&&(p=-p,m=-m))}const w=(((P=S.reduced)==null?void 0:P.keptSamples)??[]).map(T=>({x:T.x,y:T.y}));return w.length<8?null:{line:{p1:{x:f.x-p*l*.5,y:f.y-m*l*.5},p2:{x:f.x+p*l*.5,y:f.y+m*l*.5}},fitPoints:w,correctedScanlines:S.scanlines}}function Ho(n,t,e){var b;const i=n.length,r=((b=n[0])==null?void 0:b.length)??0;if(r===0||i===0)return 0;const s=Math.max(0,Math.min(r-1,t)),o=Math.max(0,Math.min(i-1,e)),a=Math.floor(s),l=Math.floor(o),u=Math.min(r-1,a+1),d=Math.min(i-1,l+1),c=s-a,g=o-l,h=n[l][a],f=n[l][u],p=n[d][a],m=n[d][u],y=h*(1-c)+f*c,x=p*(1-c)+m*c;return y*(1-g)+x*g}function Qo(n,t,e,i,r,s,o){var U;const a=i.p2.x-i.p1.x,l=i.p2.y-i.p1.y,u=Math.hypot(a,l);if(!Number.isFinite(u)||u<=1e-6)return null;const d=a/u,c=l/u,g=-c,h=d,f=(i.p1.x+i.p2.x)*.5,p=(i.p1.y+i.p2.y)*.5,m=xe(i,s+2),y=_t(Xt(m||[i.p1,i.p2],3),t,e);if(!y)return null;const x=io(n,t,e,y,0,0,(o==null?void 0:o.bayerPattern)||"RGGB",o==null?void 0:o.greenPhase,o==null?void 0:o.blackLevel),b=x.length;if((((U=x[0])==null?void 0:U.length)??0)<6||b<6)return null;const C=Math.max(8,Math.round(r*2)+1),S=Math.max(8,Math.round(s*2)+1),w=C>1?r*2/(C-1):0,P=S>1?s*2/(S-1):0,T=Array.from({length:C},()=>new Array(S).fill(0));for(let R=0;R<C;R++){const O=-r+w*R;for(let B=0;B<S;B++){const Y=-s+P*B,Z=f+O*d+Y*g,X=p+O*c+Y*h;T[R][B]=Ho(x,Z-y.x,X-y.y)}}const{gx:F,gy:k}=ro(T),_=F>=k,v=ys(T,-s,-r,P,w,_);if(v.length<8)return null;const I=v.map(R=>{const O=R.x,B=R.y;return{x:f+B*d+O*g,y:p+B*c+O*h,weight:R.weight}}),E=ve(I);if(!E)return null;let L=E.dirX,z=E.dirY;return L*d+z*c<0&&(L=-L,z=-z),{line:{p1:{x:E.pointX-L*r,y:E.pointY-z*r},p2:{x:E.pointX+L*r,y:E.pointY+z*r}},fitPoints:I.map(R=>({x:R.x,y:R.y}))}}function Pe(n){const t=Math.max(0,Math.min(1,n));return t<=.04045?t/12.92:Math.pow((t+.055)/1.055,2.4)}function er(n,t,e){if(t<=0||e<=0||n.length!==t*e)return new Uint8Array(Math.max(0,t*e));let i=1/0,r=-1/0;for(let f=0;f<n.length;f++){const p=n[f];Number.isFinite(p)&&(p<i&&(i=p),p>r&&(r=p))}if(!Number.isFinite(i)||!Number.isFinite(r)||r<=i+1e-9)return new Uint8Array(t*e);const s=1024,o=new Uint32Array(s),a=r-i;for(let f=0;f<n.length;f++){const p=Math.max(0,Math.min(1,(n[f]-i)/a)),m=Math.min(s-1,Math.max(0,Math.floor(p*(s-1))));o[m]++}const l=n.length,u=f=>{const p=l*f;let m=0;for(let y=0;y<s;y++)if(m+=o[y],m>=p)return i+y/Math.max(1,s-1)*a;return r},d=u(.01),c=u(.99),g=Math.max(1e-9,c-d),h=new Uint8Array(t*e);for(let f=0;f<n.length;f++){const p=Math.max(0,Math.min(1,(n[f]-d)/g));h[f]=Math.round(p*255)}return h}function Wn(n,t,e=0){const i=new Float32Array(n.width*n.height),r=n.data;for(let s=0,o=0;s<r.length;s+=4,o++)i[o]=_s(r,s,t,e);return er(i,n.width,n.height)}function qo(n){return Number.isFinite(n)?Math.max(0,Math.min(65535,Number(n))):0}function _s(n,t,e,i=0){let r=n[t]/255,s=n[t+1]/255,o=n[t+2]/255;e&&(r=Pe(r),s=Pe(s),o=Pe(o));const a=.2126*r+.7152*s+.0722*o;return Math.max(0,a-qo(i)/65535)}function Ko(n,t){const e=n.width,i=n.height,r=n.data;if(r.length<e*i*3)return new Uint8Array(e*i);const s=new Float32Array(e*i);for(let o=0;o<e*i;o++){const a=o*3;t!==void 0?s[o]=r[a+t]:s[o]=.2126*r[a]+.7152*r[a+1]+.0722*r[a+2]}return er(s,e,i)}function $o(n){const t=new Float32Array(n.width*n.height);for(let e=0;e<n.data.length;e++)t[e]=n.data[e];return er(t,n.width,n.height)}function Gr(n,t,e){const i=_t(t,n.width,n.height);if(!i)return null;const r=new Uint16Array(i.w*i.h*3),s=n.data;let o=0;for(let a=i.y;a<i.y+i.h;a++)for(let l=i.x;l<i.x+i.w;l++){const u=(a*n.width+l)*4;let d=s[u]/255,c=s[u+1]/255,g=s[u+2]/255;e&&(d=Pe(d),c=Pe(c),g=Pe(g)),r[o++]=Math.max(0,Math.min(65535,Math.round(d*65535))),r[o++]=Math.max(0,Math.min(65535,Math.round(c*65535))),r[o++]=Math.max(0,Math.min(65535,Math.round(g*65535)))}return{data:r,width:i.w,height:i.h}}function Yr(n,t,e,i=0){const r=_t(t,n.width,n.height);if(!r)return null;const s=new Uint16Array(r.w*r.h*3),o=n.data;let a=0;for(let l=r.y;l<r.y+r.h;l++)for(let u=r.x;u<r.x+r.w;u++){const d=(l*n.width+u)*4,c=Math.max(0,Math.min(65535,Math.round(_s(o,d,e,i)*65535)));s[a++]=c,s[a++]=c,s[a++]=c}return{data:s,width:r.w,height:r.h}}function Jo(n,t){const e=_t(t,n.width,n.height);if(!e)return null;const i=new Uint16Array(e.w*e.h),r=n.data;let s=0;for(let o=e.y;o<e.y+e.h;o++)for(let a=e.x;a<e.x+e.w;a++){const l=(o*n.width+a)*4;i[s++]=Math.max(0,Math.min(65535,Math.round((.2126*r[l]+.7152*r[l+1]+.0722*r[l+2])*257)))}return{data:i,width:e.w,height:e.h}}function Zo(n,t){const e=_t(t,n.width,n.height);if(!e)return null;const i=new Uint16Array(e.w*e.h);let r=0;for(let s=e.y;s<e.y+e.h;s++){const o=s*n.width;for(let a=e.x;a<e.x+e.w;a++)i[r++]=n.data[o+a]}return{data:i,width:e.w,height:e.h}}function Qt(n,t,e){return{x:n.x*t,y:n.y*e}}function tl(n,t,e){return{p1:Qt(n.p1,t,e),p2:Qt(n.p2,t,e)}}function Yi(n,t){const e=t(n);return{x:Number.isFinite(e.x)?e.x:n.x,y:Number.isFinite(e.y)?e.y:n.y}}function el(n,t){return n.map(e=>Yi(e,t))}function me(n,t,e,i=0,r=0){if(!n||n.length<8)return;const s=n.map(o=>({x:o.x*t-i,y:o.y*e-r})).filter(o=>Number.isFinite(o.x)&&Number.isFinite(o.y));return s.length>=8?s:void 0}function nl(n,t){return{p1:Yi(n.p1,t),p2:Yi(n.p2,t)}}function ln(n,t,e){return{p1:{x:n.p1.x-t,y:n.p1.y-e},p2:{x:n.p2.x-t,y:n.p2.y-e}}}function il(n,t,e,i){const r=Math.max(0,Math.min(n.width-1,t)),o=(Math.max(0,Math.min(n.height-1,e))*n.width+r)*4;let a=n.data[o]/255,l=n.data[o+1]/255,u=n.data[o+2]/255;return i&&(a=Pe(a),l=Pe(l),u=Pe(u)),(.2126*a+.7152*l+.0722*u)*65535}function ws(n){return n.kind==="u16-mono"}function tn(n){return n.width}function en(n){return n.height}function ei(n,t,e,i){if(ws(n)){const r=Math.max(0,Math.min(n.width-1,t)),s=Math.max(0,Math.min(n.height-1,e));return n.data[s*n.width+r]}return il(n,t,e,i)}function rl(n,t,e,i){if(ws(n)&&n.coordinateSpace==="distorted-padded"){const r=Math.round(n.paddingOffsetX??0),s=Math.round(n.paddingOffsetY??0);return ei(n,t+r,e+s,i)}return ei(n,t,e,i)}function sl(n,t,e,i,r=3){const o=[...n,{x:(n[0].x+n[1].x+n[2].x+n[3].x)*.25,y:(n[0].y+n[1].y+n[2].y+n[3].y)*.25},{x:(n[0].x+n[1].x)*.5,y:(n[0].y+n[1].y)*.5},{x:(n[1].x+n[2].x)*.5,y:(n[1].y+n[2].y)*.5},{x:(n[2].x+n[3].x)*.5,y:(n[2].y+n[3].y)*.5},{x:(n[3].x+n[0].x)*.5,y:(n[3].y+n[0].y)*.5}].map(a=>Ce(a,t)).filter(a=>Number.isFinite(a.x)&&Number.isFinite(a.y));return o.length===0?null:_t(Xt(o,r),e,i)}function nr(n,t,e,i){const r=new Map;for(let s=n.y;s<n.y+n.h;s++)for(let o=n.x;o<n.x+n.w;o++){const a=Dt({x:o,y:s},t);if(!Number.isFinite(a.x)||!Number.isFinite(a.y))continue;const l=Math.round(a.x),u=Math.round(a.y);if(l<0||u<0||l>=e||u>=i)continue;const d=r.get(u);d?(l<d.start&&(d.start=l),l>d.end&&(d.end=l)):r.set(u,{start:l,end:l})}return r}function al(n,t,e,i,r,s,o){const a=new Map,l=t.p2.x-t.p1.x,u=t.p2.y-t.p1.y,d=Math.hypot(l,u);if(!Number.isFinite(d)||d<=1e-6)return a;const c=l/d,g=u/d,h=-g,f=c,p={x:(t.p1.x+t.p2.x)*.5,y:(t.p1.y+t.p2.y)*.5},m=Math.max(1,e+1),y=Math.max(1,i+1.5);for(let x=n.y;x<n.y+n.h;x++)for(let b=n.x;b<n.x+n.w;b++){const M=b+.5,C=x+.5,S=M-p.x,w=C-p.y,P=S*c+w*g;if(!Number.isFinite(P)||Math.abs(P)>m)continue;const T=S*h+w*f;if(!Number.isFinite(T)||Math.abs(T)>y)continue;const F=Dt({x:M,y:C},r);if(!Number.isFinite(F.x)||!Number.isFinite(F.y))continue;const k=Math.round(F.x),_=Math.round(F.y);if(k<0||_<0||k>=s||_>=o)continue;const v=a.get(_);v?(k<v.start&&(v.start=k),k>v.end&&(v.end=k)):a.set(_,{start:k,end:k})}return a}function Ss(n,t,e,i,r,s){const o=new Map,a=t.p2.x-t.p1.x,l=t.p2.y-t.p1.y,u=Math.hypot(a,l);if(!Number.isFinite(u)||u<=1e-6)return o;const d=a/u,c=l/u,g=-c,h=d,f={x:(t.p1.x+t.p2.x)*.5,y:(t.p1.y+t.p2.y)*.5},p=Math.max(1,e+1),m=Math.max(1,i+1.5),y=_t(n,r,s);if(!y)return o;for(let x=y.y;x<y.y+y.h;x++)for(let b=y.x;b<y.x+y.w;b++){const M=b-f.x,C=x-f.y,S=M*d+C*c;if(!Number.isFinite(S)||Math.abs(S)>p)continue;const w=M*g+C*h;if(!Number.isFinite(w)||Math.abs(w)>m)continue;const P=o.get(x);P?(b<P.start&&(P.start=b),b>P.end&&(P.end=b)):o.set(x,{start:b,end:b})}return o}function vs(n,t,e,i){const r=new Map;for(const[s,o]of n)for(let a=o.start;a<=o.end;a++){const l=Ce({x:a,y:s},t);if(!Number.isFinite(l.x)||!Number.isFinite(l.y))continue;const u=Math.round(l.x),d=Math.round(l.y);if(u<0||d<0||u>=e||d>=i)continue;const c=r.get(d);c?(u<c.start&&(c.start=u),u>c.end&&(c.end=u)):r.set(d,{start:u,end:u})}return r}function ir(n){return Math.abs(n.k1)<1e-4&&Math.abs(n.k2)<1e-4}function ol(n){return[{x:n.x,y:n.y},{x:n.x+n.w,y:n.y},{x:n.x+n.w,y:n.y+n.h},{x:n.x,y:n.y+n.h}]}function qt(n,t,e,i,r){return Ce({x:i.x+n*t,y:i.y+n*e},r)}function rr(n,t,e,i,r){const o=qt(n,t,e,i,r),a=qt(n+1e-4,t,e,i,r);return{x:(a.x-o.x)/1e-4,y:(a.y-o.y)/1e-4}}function Ps(n,t,e,i,r,s){let o=.01;const a=d=>{const c=qt(d,t,e,i,s);return Math.hypot(c.x-r.x,c.y-r.y)},l=a(n),u=a(n+o);if(!Number.isFinite(l)||!Number.isFinite(u))return null;if(l>u){let d=n,c=n+o;for(let g=0;g<24;g++){o*=2;const h=d+o,f=a(h),p=a(c);if(!Number.isFinite(f)||!Number.isFinite(p))break;if(f>=p)return{a:d,b:h};d=c,c=h}}else{let d=n,c=n+o;for(let g=0;g<24;g++){o*=2;const h=c-o,f=a(h),p=a(d);if(!Number.isFinite(f)||!Number.isFinite(p))break;if(f>=p)return{a:h,b:c};c=d,d=h}}return{a:n-Math.max(.5,o),b:n+Math.max(.5,o)}}function ll(n,t,e=33){const i=n.p2.x-n.p1.x,r=n.p2.y-n.p1.y,s=Math.hypot(i,r);if(!Number.isFinite(s)||s<=1e-6)return[Ce(n.p1,t),Ce(n.p2,t)];const o=i/s,a=r/s,l={x:(n.p1.x+n.p2.x)*.5,y:(n.p1.y+n.p2.y)*.5},u=s*.5,d=Math.max(9,e),c=[];for(let g=0;g<d;g++){const h=d===1?.5:g/(d-1),f=-u+h*(u*2);c.push(qt(f,o,a,l,t))}return c}function cl(n,t,e,i,r,s,o=1){const a=n.p2.x-n.p1.x,l=n.p2.y-n.p1.y,u=Math.hypot(a,l);if(!Number.isFinite(u)||u<=1e-6)return null;const d=a/u,c=l/u,g={x:(n.p1.x+n.p2.x)*.5,y:(n.p1.y+n.p2.y)*.5},h=Math.max(24,Math.round(e*2)+1),f=[];for(let p=0;p<h;p++){const m=h===1?.5:p/(h-1),y=-e+m*(e*2),x=qt(y,d,c,g,t),b=rr(y,d,c,g,t),M=Math.hypot(b.x,b.y);if(!Number.isFinite(M)||M<=1e-9)continue;const C=-b.y/M,S=b.x/M;f.push({x:x.x+C*(i+o),y:x.y+S*(i+o)},{x:x.x-C*(i+o),y:x.y-S*(i+o)})}if(f.length<2){const p={p1:Ce(n.p1,t),p2:Ce(n.p2,t)},m=xe(p,i+o);return m?_t(Xt(m,o),r,s):null}return _t(Xt(f,o),r,s)}function Ze(n,t,e,i,r){return Dt({x:i.x+n*t,y:i.y+n*e},r)}function ul(n,t,e,i,r,s){let o=.01;const a=d=>{const c=Ze(d,t,e,i,s);return Math.hypot(c.x-r.x,c.y-r.y)},l=a(n),u=a(n+o);if(!Number.isFinite(l)||!Number.isFinite(u))return null;if(l>u){let d=n,c=n+o;for(let g=0;g<24;g++){o*=2;const h=d+o,f=a(h),p=a(c);if(!Number.isFinite(f)||!Number.isFinite(p))break;if(f>=p)return{a:d,b:h};d=c,c=h}}else{let d=n,c=n+o;for(let g=0;g<24;g++){o*=2;const h=c-o,f=a(h),p=a(d);if(!Number.isFinite(f)||!Number.isFinite(p))break;if(f>=p)return{a:h,b:c};c=d,d=h}}return{a:n-Math.max(.5,o),b:n+Math.max(.5,o)}}function sr(n,t,e){const i=(t.x-n.x)*(t.y-e.y)-(t.x-e.x)*(t.y-n.y);if(Math.abs(i)<=1e-12)return .5*(n.x+e.x);const r=(t.x-n.x)*(t.x-n.x)*(t.y-e.y)-(t.x-e.x)*(t.x-e.x)*(t.y-n.y),s=t.x-.5*r/i;return Number.isFinite(s)?s:.5*(n.x+e.x)}function hl(n,t,e,i,r){const o=Ze(n,t,e,i,r),a=Ze(n+1e-4,t,e,i,r);return{x:(a.x-o.x)/1e-4,y:(a.y-o.y)/1e-4}}function dl(n,t,e,i,r,s,o=!1,a){const l=[],u=[],d=a?Et*2:Et,c=Math.max(1,Math.min(s,d)),g=i.p2.x-i.p1.x,h=i.p2.y-i.p1.y,f=Math.hypot(g,h);if(!Number.isFinite(f)||f<=1e-6)return null;const p=g/f,m=h/f,y={x:(i.p1.x+i.p2.x)*.5,y:(i.p1.y+i.p2.y)*.5},x={p1:Dt(i.p1,e),p2:Dt(i.p2,e)},b=x.p2.x-x.p1.x,M=x.p2.y-x.p1.y,C=Math.hypot(b,M);if(!Number.isFinite(C)||C<=1e-6)return null;const S=b/C,w=M/C,P=-w,T=S,F={x:(x.p1.x+x.p2.x)*.5,y:(x.p1.y+x.p2.y)*.5},k=_t(a||Xt(xe(i,s+2)??[i.p1,i.p2],2),n.width,n.height);if(!k)return null;const _=al(k,i,r,c,e,tn(t),en(t));if(_.size===0)return null;const v=!ir(e);for(const[R,O]of _)if(!(R<0||R>=en(t)))for(let B=O.start;B<=O.end;B++){if(B<0||B>=tn(t))continue;const Y={x:B,y:R};let Z,X;if(v){const J=Ce(Y,e);if(!Number.isFinite(J.x)||!Number.isFinite(J.y)||Math.round(J.x)<0||Math.round(J.x)>=n.width||Math.round(J.y)<0||Math.round(J.y)>=n.height)continue;const W=J.x-y.x,tt=J.y-y.y,lt=W*p+tt*m;if(!Number.isFinite(lt))continue;Z=lt,X=W*-m+tt*p;const D=ul(lt,p,m,y,Y,e);if(!D)continue;const $=.5*(D.a+D.b),K=Ze(D.a,p,m,y,e),A=Ze($,p,m,y,e),Q=Ze(D.b,p,m,y,e),G=sr({x:D.a,y:Math.hypot(K.x-Y.x,K.y-Y.y)},{x:$,y:Math.hypot(A.x-Y.x,A.y-Y.y)},{x:D.b,y:Math.hypot(Q.x-Y.x,Q.y-Y.y)});if(!Number.isFinite(G))continue;Z=G;const V=hl(G,p,m,y,e),rt=Math.hypot(V.x,V.y);if(!Number.isFinite(rt)||rt<=1e-9)continue;const st=V.x/rt,dt=-(V.y/rt),ft=st,et=Ze(G,p,m,y,e);X=(Y.x-et.x)*dt+(Y.y-et.y)*ft}else{const J=Y.x-F.x,W=Y.y-F.y;Z=J*S+W*w,X=J*P+W*T}!Number.isFinite(Z)||Math.abs(Z)>r||!Number.isFinite(X)||Math.abs(X)>c||(l.push(X),u.push(ei(t,B,R,o)))}if(l.length<8)return null;const I=Dt(i.p1,e),E=Dt(i.p2,e),L=E.x-I.x,z=E.y-I.y,U=Math.abs(L)>=Math.abs(z)?1:2;return yn(l,u,U,d)}function fl(n,t,e,i,r,s,o=!1,a,l,u,d=!1){const c=[],g=[],h=a?Et*2:Et,f=Math.max(1,Math.min(s,h)),p=i.p2.x-i.p1.x,m=i.p2.y-i.p1.y,y=Math.hypot(p,m);if(!Number.isFinite(y)||y<=1e-6)return null;const x=p/y,b=m/y,M=-b,C=x,S={x:(i.p1.x+i.p2.x)*.5,y:(i.p1.y+i.p2.y)*.5},w=_t(a||Xt(xe(i,h*4+2)??[i.p1,i.p2],2),n.width,n.height);if(!w)return null;const P=l??(u?nr(_t(u,tn(t),en(t))??u,e,n.width,n.height):Ss(w,i,Math.max(1,r),f*4+.5,n.width,n.height));if(P.size===0)return null;const T=vs(P,e,tn(t),en(t));if(T.size===0)return null;const F=!ir(e);for(const[_,v]of T)for(let I=v.start;I<=v.end;I++){const E={x:I,y:_},L=Dt(E,e);if(!Number.isFinite(L.x)||!Number.isFinite(L.y)||Math.round(L.x)<0||Math.round(L.x)>=n.width||Math.round(L.y)<0||Math.round(L.y)>=n.height)continue;const z=L.x-S.x,U=L.y-S.y,R=z*x+U*b;let O=z*M+U*C;if(F){const B=Ps(R,x,b,S,E,e);if(!B)continue;const Y=.5*(B.a+B.b),Z=qt(B.a,x,b,S,e),X=qt(Y,x,b,S,e),J=qt(B.b,x,b,S,e),W=sr({x:B.a,y:Math.hypot(Z.x-E.x,Z.y-E.y)},{x:Y,y:Math.hypot(X.x-E.x,X.y-E.y)},{x:B.b,y:Math.hypot(J.x-E.x,J.y-E.y)});if(!Number.isFinite(W))continue;const tt=rr(W,x,b,S,e),lt=Math.hypot(tt.x,tt.y);if(!Number.isFinite(lt)||lt<=1e-9)continue;const D=tt.x/lt,K=-(tt.y/lt),A=D,Q=qt(W,x,b,S,e);O=(E.x-Q.x)*K+(E.y-Q.y)*A}!Number.isFinite(R)||Math.abs(R)>Math.max(1,r)||!Number.isFinite(O)||Math.abs(O)>f||(c.push(O),g.push(rl(t,I,_,o)))}if(c.length<8)return null;const k=Math.abs(p)>=Math.abs(m)?1:2;return d?kn(c,g,k,h):yn(c,g,k,h)}function pl(n,t,e,i,r,s){const o=n.width,a=n.height,l=(n.bayerPattern||"RGGB").toUpperCase(),u=s!=null&&s.correctedRect?Et*2:Et,d=Math.max(1,Math.min(r,u)),c=(s==null?void 0:s.restrictToStrip)??!0,g=e.p2.x-e.p1.x,h=e.p2.y-e.p1.y,f=Math.hypot(g,h);if(!Number.isFinite(f)||f<=1e-6)return null;const p=g/f,m=h/f,y=-m,x=p,b={x:(e.p1.x+e.p2.x)*.5,y:(e.p1.y+e.p2.y)*.5},M={p1:{x:b.x-p*Math.max(1,i),y:b.y-m*Math.max(1,i)},p2:{x:b.x+p*Math.max(1,i),y:b.y+m*Math.max(1,i)}},C=xe(M,d+2),S=(s!=null&&s.fixedRawRect?_t(s.fixedRawRect,o,a):null)??(s!=null&&s.correctedRect?Zi(s.correctedRect,t,o,a):null)??(C?sl(C,t,o,a,2):null);if(!S)return null;const w=[],P=[];for(let k=S.y;k<S.y+S.h;k++){const _=k*o;for(let v=S.x;v<S.x+S.w;v++){if(!vt(v,k,l,s==null?void 0:s.greenPhase))continue;const I=Dt({x:v,y:k},t);if(!Number.isFinite(I.x)||!Number.isFinite(I.y))continue;const E=I.x-b.x,L=I.y-b.y,z=E*p+L*m;if(!Number.isFinite(z)||c&&Math.abs(z)>Math.max(1,i))continue;const U=E*y+L*x;if(!Number.isFinite(U)||c&&Math.abs(U)>d)continue;w.push(U);let R;R=Math.max(0,n.data[_+v]-gn(s==null?void 0:s.blackLevel,v,k)),P.push(R)}}if(w.length<8)return null;const T=Math.abs(g)>=Math.abs(h)?1:2,F=Math.max(2,(s==null?void 0:s.shortSidePxOverride)??(c?d*2:Math.min(S.w,S.h)));return ci(w,P,F,s==null?void 0:s.manualBinSize,T,s==null?void 0:s.preferAutoPerEdgeBin,!1,!!(s!=null&&s.forceLegacyModel))}function Cs(n,t,e,i,r,s=!1,o,a,l,u=!1){if(a&&l){const P=l.p2.x-l.p1.x,T=l.p2.y-l.p1.y,F=Math.hypot(P,T);if(Number.isFinite(F)&&F>1e-6)return dl(a,n,t,l,Math.max(1,F*.5),r,s,o)}const d=[],c=[],g=Et,h=e.p2.x-e.p1.x,f=e.p2.y-e.p1.y,p=Math.hypot(h,f);if(!Number.isFinite(p)||p<=1e-6)return null;const m=h/p,y=f/p,x={x:(e.p1.x+e.p2.x)*.5,y:(e.p1.y+e.p2.y)*.5},b=(o?_t(o,tn(n),en(n)):null)??cl(e,t,i+1,r+1,tn(n),en(n),1);if(!b)return null;const M=nr(b,t,tn(n),en(n));if(M.size===0)return null;const C=-y,S=m;for(const[P,T]of M)for(let F=T.start;F<=T.end;F++){const k={x:F,y:P};let _=(k.x-x.x)*m+(k.y-x.y)*y,v=(k.x-x.x)*C+(k.y-x.y)*S;!Number.isFinite(_)||Math.abs(_)>i+1||!Number.isFinite(v)||Math.abs(v)>=g||(d.push(v),c.push(ei(n,F,P,s)))}if(d.length<8)return null;const w=Math.abs(h)>=Math.abs(f)?1:2;return u?kn(d,c,w,g):yn(d,c,w,g)}function ml(n,t,e){var g;const i=e.sourceMode??(t.isThreePlane?"three-plane":"rggb-raw"),r=e.useQuadraticProjection!==!1,s=!!e.forceRenderedMeasurement,o=n.width,a=n.height,l=e.threePlaneChannel,u=tr(e.detectionTuning),d=e.monochromeBlackLevel??0;if(i==="rggb-raw"&&!s){if(!t||t.isThreePlane)return null;const h=So(t,e.greenPhase),f=o/Math.max(1,t.width),p=a/Math.max(1,t.height);return{sourceMode:i,detectionGray:h,detectionWidth:t.width,detectionHeight:t.height,detectToMeasureX:1,detectToMeasureY:1,detectToDisplayX:f,detectToDisplayY:p,measureToDisplayX:f,measureToDisplayY:p,detectPointToDisplay:m=>Qt(m,f,p),measurePointToDisplay:m=>Qt(m,f,p),displayPointToDetect:m=>Qt(m,1/Math.max(1e-9,f),1/Math.max(1e-9,p)),measureUsesDisplayLine:!1,measureWidth:t.width,measureHeight:t.height,refineLine:(m,y,x)=>Qo(t.data,t.width,t.height,{p1:m,p2:y},x*.5,Math.max(4,x*.2),{greenPhase:e.greenPhase,bayerPattern:t.bayerPattern})||ge(h,t.width,t.height,m,y,x*.5,Math.max(4,x*.2)),measureEdge:(m,y,x,b,M)=>Se(t.data,t.width,t.height,m,y,x,b,{greenOnly:!0,greenPhase:e.greenPhase,bayerPattern:t.bayerPattern,blackLevel:e.blackLevel??void 0,preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?me(M==null?void 0:M.fitPoints,1,1):void 0})}}if(s){const h=!!e.distortionCurveApplied&&!!e.distortionModel,f=i==="rggb-raw"&&!!e.distortionCorrected&&!!e.distortionModel&&!t.isThreePlane,p=!!e.distortionCorrected&&!!e.distortionModel&&!!e.distortionOriginalSamplingPlane,m=!!e.distortionCorrected&&!!e.distortionSamplingPlane,y=n,x=Wn(y,!!e.sfrHasGamma,i==="unmix-bw"?d:0);return{sourceMode:i,detectionGray:x,detectionWidth:y.width,detectionHeight:y.height,detectToMeasureX:1,detectToMeasureY:1,detectToDisplayX:1,detectToDisplayY:1,measureToDisplayX:1,measureToDisplayY:1,detectPointToDisplay:b=>b,measurePointToDisplay:b=>b,displayPointToDetect:b=>b,measureUsesDisplayLine:!1,measureWidth:y.width,measureHeight:y.height,refineLine:(b,M,C)=>(p?Wo(x,n.width,n.height,b,M,Et):null)||ge(x,n.width,n.height,b,M,C*.5,Math.max(4,C*.2)),measureEdge:(b,M,C,S,w)=>{const P=e.distortionModel?Zi(b,e.distortionModel,t.width,t.height):null;if(f){const _={p1:Dt(M.p1,e.distortionModel),p2:Dt(M.p2,e.distortionModel)},v=Math.hypot(_.p2.x-_.p1.x,_.p2.y-_.p1.y),I=Math.max(2,v*.5*u.sampleHalfWidthRatio);return Xa(t,e.distortionModel,_,Math.max(1,v*.5),I,{greenPhase:e.greenPhase,blackLevel:e.blackLevel??void 0,correctedRect:b})}if(f)return pl(t,e.distortionModel,M,C,S,{greenPhase:e.greenPhase,blackLevel:e.blackLevel??void 0,correctedRect:b,fixedRawRect:P,preferAutoPerEdgeBin:!0});if(p)return fl(n,e.distortionOriginalSamplingPlane,e.distortionModel,M,C,S,!!e.sfrHasGamma,b,(w==null?void 0:w.correctedScanlines)??null,P);if(h){const _={p1:Dt(M.p1,e.distortionModel),p2:Dt(M.p2,e.distortionModel)},v=Math.hypot(_.p2.x-_.p1.x,_.p2.y-_.p1.y);return Cs(e.distortionOriginalSamplingPlane??e.distortionSamplingPlane??e.distortionSamplingImage??n,e.distortionModel,_,Math.max(1,v*.5),S,!!e.sfrHasGamma,b,e.distortionBaseImage??n,M)}if(m){const _=Zo(e.distortionSamplingPlane,b);if(!_)return null;const v=ln(M,b.x,b.y);return Se(_.data,_.width,_.height,{x:0,y:0,w:_.width,h:_.height},v,C,S,{preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?me(w==null?void 0:w.fitPoints,1,1,b.x,b.y):void 0})}const F=i==="unmix-bw"?Yr(n,b,!!e.sfrHasGamma,d):Gr(n,b,!!e.sfrHasGamma);if(!F)return null;const k=ln(M,b.x,b.y);return Se(F.data,F.width,F.height,{x:0,y:0,w:F.width,h:F.height},k,C,S,{isThreePlane:!0,threePlaneChannel:void 0,preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?me(w==null?void 0:w.fitPoints,1,1,b.x,b.y):void 0})}}}if(i==="three-plane"){if(t.isThreePlane&&!e.sfrHasGamma){const f=Ko(t,l),p=o/Math.max(1,t.width),m=a/Math.max(1,t.height);return{sourceMode:i,detectionGray:f,detectionWidth:t.width,detectionHeight:t.height,detectToMeasureX:1,detectToMeasureY:1,detectToDisplayX:p,detectToDisplayY:m,measureToDisplayX:p,measureToDisplayY:m,detectPointToDisplay:y=>Qt(y,p,m),measurePointToDisplay:y=>Qt(y,p,m),displayPointToDetect:y=>Qt(y,1/Math.max(1e-9,p),1/Math.max(1e-9,m)),measureUsesDisplayLine:!1,measureWidth:t.width,measureHeight:t.height,refineLine:(y,x,b)=>ge(f,t.width,t.height,y,x,b*.5,Math.max(4,b*.2)),measureEdge:(y,x,b,M,C)=>Se(t.data,t.width,t.height,y,x,b,M,{isThreePlane:!0,threePlaneChannel:l,preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?me(C==null?void 0:C.fitPoints,1,1):void 0})}}const h=Wn(n,!!e.sfrHasGamma);return{sourceMode:i,detectionGray:h,detectionWidth:n.width,detectionHeight:n.height,detectToMeasureX:1,detectToMeasureY:1,detectToDisplayX:1,detectToDisplayY:1,measureToDisplayX:1,measureToDisplayY:1,detectPointToDisplay:f=>f,measurePointToDisplay:f=>f,displayPointToDetect:f=>f,measureUsesDisplayLine:!1,measureWidth:n.width,measureHeight:n.height,refineLine:(f,p,m)=>ge(h,n.width,n.height,f,p,m*.5,Math.max(4,m*.2)),measureEdge:(f,p,m,y,x)=>{const b=Gr(n,f,!!e.sfrHasGamma);if(!b)return null;const M=ln(p,f.x,f.y);return Se(b.data,b.width,b.height,{x:0,y:0,w:b.width,h:b.height},M,m,y,{isThreePlane:!0,threePlaneChannel:l,preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?me(x==null?void 0:x.fitPoints,1,1,f.x,f.y):void 0})}}}if(i==="unmix-bw"){if(t&&!t.isThreePlane&&e.displaySettings){const f=wa(t,e.displaySettings,e.blackLevel??e.monochromeBlackLevel??void 0);if(f){const p=$o(f),m=o/Math.max(1,t.width),y=a/Math.max(1,t.height);return{sourceMode:i,detectionGray:p,detectionWidth:t.width,detectionHeight:t.height,detectToMeasureX:1,detectToMeasureY:1,detectToDisplayX:m,detectToDisplayY:y,measureToDisplayX:m,measureToDisplayY:y,detectPointToDisplay:x=>Qt(x,m,y),measurePointToDisplay:x=>Qt(x,m,y),displayPointToDetect:x=>Qt(x,1/Math.max(1e-9,m),1/Math.max(1e-9,y)),measureUsesDisplayLine:!1,measureWidth:t.width,measureHeight:t.height,refineLine:(x,b,M)=>ge(p,t.width,t.height,x,b,M*.5,Math.max(4,M*.2)),measureEdge:(x,b,M,C,S)=>Se(f.data,f.width,f.height,x,b,M,C,{preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?me(S==null?void 0:S.fitPoints,1,1):void 0})}}}const h=Wn(n,!!e.sfrHasGamma,d);return{sourceMode:i,detectionGray:h,detectionWidth:n.width,detectionHeight:n.height,detectToMeasureX:1,detectToMeasureY:1,detectToDisplayX:1,detectToDisplayY:1,measureToDisplayX:1,measureToDisplayY:1,detectPointToDisplay:f=>f,measurePointToDisplay:f=>f,displayPointToDetect:f=>f,measureUsesDisplayLine:!1,measureWidth:n.width,measureHeight:n.height,refineLine:(f,p,m)=>ge(h,n.width,n.height,f,p,m*.5,Math.max(4,m*.2)),measureEdge:(f,p,m,y,x)=>{const b=Yr(n,f,!!e.sfrHasGamma,d);if(!b)return null;const M=ln(p,f.x,f.y);return Se(b.data,b.width,b.height,{x:0,y:0,w:b.width,h:b.height},M,m,y,{isThreePlane:!0,threePlaneChannel:void 0,preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?me(x==null?void 0:x.fitPoints,1,1,f.x,f.y):void 0})}}}const c=Wn(n,!1);if(t&&!t.isThreePlane&&((g=e.displaySettings)==null?void 0:g.renderMode)==="advanced-zero-dep"&&e.displaySettings.advancedZeroDep){const h=o/Math.max(1,t.width),f=a/Math.max(1,t.height);return{sourceMode:i,detectionGray:c,detectionWidth:n.width,detectionHeight:n.height,detectToMeasureX:t.width/Math.max(1,n.width),detectToMeasureY:t.height/Math.max(1,n.height),detectToDisplayX:1,detectToDisplayY:1,measureToDisplayX:h,measureToDisplayY:f,detectPointToDisplay:p=>p,measurePointToDisplay:p=>Qt(p,h,f),displayPointToDetect:p=>p,measureUsesDisplayLine:!1,measureWidth:t.width,measureHeight:t.height,refineLine:(p,m,y)=>ge(c,n.width,n.height,p,m,y*.5,Math.max(4,y*.2)),measureEdge:(p,m,y,x,b)=>{const M=_a(t,p,e.displaySettings);if(!M||M.width<8||M.height<8)return null;const C=ln(m,p.x,p.y);return Se(M.data,M.width,M.height,{x:0,y:0,w:M.width,h:M.height},C,y,x,{preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?me(b==null?void 0:b.fitPoints,t.width/Math.max(1,n.width),t.height/Math.max(1,n.height),p.x,p.y):void 0})}}}if(t&&!t.isThreePlane){const h=o/Math.max(1,t.width),f=a/Math.max(1,t.height);return{sourceMode:i,detectionGray:c,detectionWidth:n.width,detectionHeight:n.height,detectToMeasureX:t.width/Math.max(1,n.width),detectToMeasureY:t.height/Math.max(1,n.height),detectToDisplayX:1,detectToDisplayY:1,measureToDisplayX:h,measureToDisplayY:f,detectPointToDisplay:p=>p,measurePointToDisplay:p=>Qt(p,h,f),displayPointToDetect:p=>p,measureUsesDisplayLine:!1,measureWidth:t.width,measureHeight:t.height,refineLine:(p,m,y)=>ge(c,n.width,n.height,p,m,y*.5,Math.max(4,y*.2)),measureEdge:(p,m,y,x,b)=>Se(t.data,t.width,t.height,p,m,y,x,{blackLevel:e.blackLevel??void 0,preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?me(b==null?void 0:b.fitPoints,t.width/Math.max(1,n.width),t.height/Math.max(1,n.height)):void 0})}}return{sourceMode:i,detectionGray:c,detectionWidth:n.width,detectionHeight:n.height,detectToMeasureX:1,detectToMeasureY:1,detectToDisplayX:1,detectToDisplayY:1,measureToDisplayX:1,measureToDisplayY:1,detectPointToDisplay:h=>h,measurePointToDisplay:h=>h,displayPointToDetect:h=>h,measureUsesDisplayLine:!1,measureWidth:n.width,measureHeight:n.height,refineLine:(h,f,p)=>ge(c,n.width,n.height,h,f,p*.5,Math.max(4,p*.2)),measureEdge:(h,f,p,m,y)=>{const x=Jo(n,h);if(!x)return null;const b=ln(f,h.x,h.y);return mo(x.data,x.width,x.height,{x:0,y:0,w:x.width,h:x.height},b,p,m,{blackLevel:e.monochromeBlackLevel??void 0,preferAutoPerEdgeBin:!0,disableQuadraticProjection:!r,quadraticFitPoints:r?me(y==null?void 0:y.fitPoints,1,1,h.x,h.y):void 0})}}}function gl(n,t,e){var u,d,c,g,h,f,p,m;if(!n||!t)return[];(u=e.onProgress)==null||u.call(e,"Preparing source context...",0);const i=ml(n,t,e);if(!i)return[];(d=e.onProgress)==null||d.call(e,"Preparing source context...",.08);const r=tr(e.detectionTuning),s=Math.min(1e3,Math.max(1,e.maxRegions??1e3)),o=Math.max(4,e.maxEdges??s*4);(c=e.onProgress)==null||c.call(e,"Detecting candidates...",.12);const a=Bo(i.detectionGray,i.detectionWidth,i.detectionHeight,s,e.detectionTuning,(y,x)=>{var b;(b=e.onProgress)==null||b.call(e,y,.12+.08*Math.max(0,Math.min(1,x)))});(g=e.onProgress)==null||g.call(e,"Detecting candidates...",.2);const l=[];for(let y=0;y<a.length;y++){const x=a.length<=0?1:y/a.length;if((h=e.onProgress)==null||h.call(e,`Measuring edges: region ${y+1}/${a.length}`,.2+.72*Math.min(1,x)),l.length>=o)break;const b=a[y],M=b.corners,C=`auto-region-${y+1}`;for(let S=0;S<4&&((f=e.onProgress)==null||f.call(e,`Measuring edges: region ${y+1}/${a.length}, edge ${S+1}/4`,.2+.72*Math.min(1,(y+S/4)/Math.max(1,a.length))),!(l.length>=o));S=S+1){const w=M[S],P=M[(S+1)%4],T=P.x-w.x,F=P.y-w.y,k=Math.hypot(T,F);if(!Number.isFinite(k)||k<24)continue;const _=.125,v={x:w.x+T*_,y:w.y+F*_},I={x:P.x-T*_,y:P.y-F*_},E=Math.hypot(I.x-v.x,I.y-v.y);if(!Number.isFinite(E)||E<12)continue;const L=i.refineLine(v,I,E),z=(L!=null&&L.fitPoints?Rr(L.fitPoints):null)||(L==null?void 0:L.line)||{p1:v,p2:I},U=tl(z,i.detectToMeasureX,i.detectToMeasureY),R=el((L==null?void 0:L.fitPoints)??[],i.detectPointToDisplay),O=(R.length>=2?Rr(R):null)||nl(U,i.measurePointToDisplay),B=i.measureUsesDisplayLine?O:U,Y=B.p2.x-B.p1.x,Z=B.p2.y-B.p1.y,X=Math.hypot(Y,Z);if(!Number.isFinite(X)||X<=1e-6)continue;const J=O.p2.x-O.p1.x,W=O.p2.y-O.p1.y,tt=Math.hypot(J,W);if(!Number.isFinite(tt)||tt<=1e-6)continue;const lt=!!e.distortionCurveApplied&&!!e.distortionModel,D=J/tt,$=W/tt;let K=$,A=-D;const Q=(O.p1.x+O.p2.x)*.5,G=(O.p1.y+O.p2.y)*.5,V=i.detectPointToDisplay({x:b.centerX,y:b.centerY}),rt=V.x,st=V.y;(Q-rt)*K+(G-st)*A<0&&(K=-K,A=-A);const it=X*.5,dt=Math.max(2,X*r.sampleHalfWidthRatio),ft=Math.max(2,tt*r.sampleHalfWidthRatio),et=lt?{p1:Dt(O.p1,e.distortionModel),p2:Dt(O.p2,e.distortionModel)}:void 0,Bt=et?Math.max(1,Math.hypot(et.p2.x-et.p1.x,et.p2.y-et.p1.y)*.5):it,Gt=lt?O:U,j=lt?ft:dt,ct=xe(Gt,j);if(!ct)continue;const at=_t(Xt(ct,2),lt?n.width:i.measureWidth,lt?n.height:i.measureHeight);if(!at)continue;const pt=e.distortionCorrected&&e.distortionModel&&i.sourceMode==="rggb-raw"?Zi(at,e.distortionModel,t.width,t.height):null,mt=lt?Cs(e.distortionSamplingPlane??e.distortionSamplingImage??n,e.distortionModel,et,Bt,j,!!e.sfrHasGamma,at,e.distortionBaseImage??n,Gt):i.measureEdge(at,B,it,j,L);if(!mt||(mt.autoLikeUsed=!0,!fo(mt,e.useDeshading,0)))continue;const te=e.useNR?-1:12,Ut=Mo([mt],te,null,e.useDeshading,0,!0);if(!Ut||Ut.mtf50===null||!po(Ut.lsfCropped))continue;const Mt=et?ll(et,e.distortionModel,Math.max(21,Math.round(tt*.5))):mt.quadraticProjectionUsed?ho(R,O,Math.max(21,Math.round(tt*.5))):void 0,le=et&&Mt&&Mt.length>=2?Xt(Mt,ft+2):null,Kt=le?ol(le):xe(O,ft);if(!Kt)continue;const ke=e.distortionCorrected?at:le??Xt(Kt,2);let q={x:Q+K*(ft+12),y:G+A*(ft+12)},Ft=Er(D,$);if(Mt&&Mt.length>=3){const At=Math.floor(Mt.length/2),xt=Mt[Math.max(0,At-1)],wt=Mt[Math.min(Mt.length-1,At+1)],ee=Mt[At],$t=wt.x-xt.x,xn=wt.y-xt.y,Be=Math.hypot($t,xn);if(Be>1e-6){const Ue=xn/Be,Oe=-$t/Be;Ft=Er($t/Be,xn/Be);const Xe={x:ee.x-rt,y:ee.y-st},Ge=Xe.x*Ue+Xe.y*Oe>=0?1:-1;q={x:ee.x+Ue*Ge*(ft+12),y:ee.y+Oe*Ge*(ft+12)}}}l.push({id:`${C}-edge-${S+1}`,regionId:C,sourceMode:i.sourceMode,edgeIndex:S,label:Ut.mtf50.toFixed(3),mtf50:Ut.mtf50,angle:Ft,orientation:mt.orientation,edgeData:mt,sourceRect:ke,rawSourceRect:(i.sourceMode==="rggb-raw"?pt??at:pt)??void 0,quad:Kt,line:O,originalLine:U,curveBaseLine:et,curvePoints:Mt,labelPoint:q,ridgePoints:R,outerSideMeans:b.outerSideMeans,outerSideQuads:b.outerSideQuads,distortionCorrected:e.distortionCorrected??!1})}}return(p=e.onProgress)==null||p.call(e,"Finalizing results...",.98),(m=e.onProgress)==null||m.call(e,"Finalizing results...",1),l}const yl=n=>!n.blackLevels||n.blackLevels.length<4?null:[Number(n.blackLevels[0])||0,Number(n.blackLevels[1])||0,Number(n.blackLevels[2])||0,Number(n.blackLevels[3])||0],Ei=(n,t)=>{t instanceof ArrayBuffer&&(n.includes(t)||n.push(t))};self.onmessage=async n=>{var o,a;const{id:t,buffer:e,detect:i,options:r}=n.data,s=performance.now();try{const l=performance.now(),u=await Ma(e),d=performance.now()-l;let c=0,g=[];if(i&&!u.isXTrans){const f=u.isThreePlane?"three-plane":"rggb-raw",p=performance.now();g=gl({width:u.width,height:u.height},u,{...r,sourceMode:f,forceRenderedMeasurement:!1,blackLevel:(r==null?void 0:r.blackLevel)??yl(u),onProgress:(m,y)=>{self.postMessage({id:t,type:"progress",stage:m,progress:y})}}),c=performance.now()-p}const h=[];Ei(h,e),Ei(h,(o=u.data)==null?void 0:o.buffer),Ei(h,(a=u.floatData)==null?void 0:a.buffer),self.postMessage({id:t,type:"result",success:!0,raw:u,rawFileBuffer:e,measurements:g,timings:{decodeMs:d,detectMs:c,totalMs:performance.now()-s}},h)}catch(l){self.postMessage({id:t,type:"result",success:!1,error:(l==null?void 0:l.message)||String(l)})}};
