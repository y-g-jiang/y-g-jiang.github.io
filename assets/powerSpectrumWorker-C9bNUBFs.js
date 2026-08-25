var Sn=Object.defineProperty;var kn=(r,e,t)=>e in r?Sn(r,e,{enumerable:!0,configurable:!0,writable:!0,value:t}):r[e]=t;var z=(r,e,t)=>kn(r,typeof e!="symbol"?e+"":e,t);class He{constructor(e){z(this,"size");z(this,"isPowerOfTwo");z(this,"_real");z(this,"_imag");z(this,"_scratch",null);z(this,"_rev",null);z(this,"_m",0);z(this,"_internalFFT",null);z(this,"_chirpReal",null);z(this,"_chirpImag",null);z(this,"_bReal",null);z(this,"_bImag",null);z(this,"_hanning",null);z(this,"_windowSumSq",0);this.size=e,this.isPowerOfTwo=(e&e-1)===0&&e>0,this._real=new Float32Array(e),this._imag=new Float32Array(e),this.isPowerOfTwo?this.initRadix2():this.initBluestein()}initRadix2(){const e=this.size,t=Math.log2(e);this._rev=new Uint32Array(e);for(let n=0;n<e;n++){let i=0,s=n;for(let a=0;a<t;a++)i=i<<1|s&1,s>>>=1;this._rev[n]=i}}initBluestein(){const e=this.size;this._m=Math.pow(2,Math.ceil(Math.log2(2*e-1))),this._internalFFT=new He(this._m),this._chirpReal=new Float32Array(e),this._chirpImag=new Float32Array(e);for(let i=0;i<e;i++){const s=-Math.PI*(i*i)/e;this._chirpReal[i]=Math.cos(s),this._chirpImag[i]=Math.sin(s)}const t=new Float32Array(this._m),n=new Float32Array(this._m);for(let i=0;i<e;i++)t[i]=this._chirpReal[i],n[i]=-this._chirpImag[i];for(let i=1;i<e;i++)t[this._m-i]=t[i],n[this._m-i]=n[i];this._internalFFT.transform(t,n),this._bReal=new Float32Array(this._internalFFT._real),this._bImag=new Float32Array(this._internalFFT._imag)}initHanning(){if(this._hanning)return;const e=this.size;this._hanning=new Float32Array(e);let t=0;for(let n=0;n<e;n++){const i=.5*(1-Math.cos(2*Math.PI*n/(e-1)));this._hanning[n]=i,t+=i*i}this._windowSumSq=t}transform(e,t){this.isPowerOfTwo?this.transformRadix2(e,t):this.transformBluestein(e,t)}transformRadix2(e,t){const n=this.size,i=this._rev,s=this._real,a=this._imag;if(e===s)for(let o=0;o<n;o++){const l=i[o];if(o<l){const h=s[o],c=a[o];s[o]=s[l],a[o]=a[l],s[l]=h,a[l]=c}}else for(let o=0;o<n;o++){const l=i[o];s[o]=e[l],a[o]=t?t[l]:0}for(let o=2;o<=n;o*=2){const l=o/2,h=-2*Math.PI/o,c=Math.cos(h),d=Math.sin(h);for(let u=0;u<n;u+=o){let p=1,_=0;for(let g=0;g<l;g++){const y=u+g,b=u+g+l,v=p*s[b]-_*a[b],A=p*a[b]+_*s[b],T=s[y],P=a[y];s[y]=T+v,a[y]=P+A,s[b]=T-v,a[b]=P-A;const O=p*c-_*d,E=p*d+_*c;p=O,_=E}}}}transformBluestein(e,t){const n=this.size,i=this._m,s=this._internalFFT,a=s._real,o=s._imag;a.fill(0),o.fill(0);for(let d=0;d<n;d++){const u=e[d],p=t?t[d]:0,_=this._chirpReal[d],g=this._chirpImag[d];a[d]=u*_-p*g,o[d]=u*g+p*_}s.transformRadix2(a,o);for(let d=0;d<i;d++){const u=s._real[d],p=s._imag[d],_=this._bReal[d],g=this._bImag[d];s._real[d]=u*_-p*g,s._imag[d]=u*g+p*_}const l=s._real,h=s._imag;for(let d=0;d<i;d++)h[d]=-h[d];s.transformRadix2(l,h);const c=1/i;for(let d=0;d<n;d++){const u=s._real[d]*c,p=-s._imag[d]*c,_=this._chirpReal[d],g=this._chirpImag[d];this._real[d]=u*_-p*g,this._imag[d]=u*g+p*_}}calculateSpectrum(e,t,n=!1){const i=this.size;let s=0;for(let c=0;c<i;c++)s+=e[c];const a=s/i;this._scratch||(this._scratch=new Float32Array(i));const o=this._scratch;if(n){this.initHanning();const c=this._hanning;for(let d=0;d<i;d++)o[d]=(e[d]-a)*c[d]}else for(let c=0;c<i;c++)o[c]=e[c]-a;this.transform(o);const l=t.length;let h=1/i;n&&this._windowSumSq>0&&(h=1/this._windowSumSq);for(let c=0;c<l;c++){const d=this._real[c],u=this._imag[c];t[c]+=(d*d+u*u)*h}}calculateSpectrumWindow(e,t,n,i=!1){const s=this.size;let a=0;for(let d=0;d<s;d++)a+=e[t+d];const o=a/s;this._scratch||(this._scratch=new Float32Array(s));const l=this._scratch;if(i){this.initHanning();const d=this._hanning;for(let u=0;u<s;u++)l[u]=(e[t+u]-o)*d[u]}else for(let d=0;d<s;d++)l[d]=e[t+d]-o;this.transform(l);const h=n.length;let c=1/s;i&&this._windowSumSq>0&&(c=1/this._windowSumSq);for(let d=0;d<h;d++){const u=this._real[d],p=this._imag[d];n[d]+=(u*u+p*p)*c}}}const Cn={"Sony ILCE-7RM5":"0.82 -0.2976 -0.0719 -0.4296 1.2053 0.2532 -0.0429 0.1282 0.5774"};let at=null;async function vn(r){return at||(at=(async()=>{if(typeof window.loadPyodide!="function")throw new Error("Pyodide missing: window.loadPyodide not found.");const e=await window.loadPyodide();return await e.loadPackage("numpy"),e})()),at}let ot=null,Rt=!1;async function bt(){var t;ot||(ot=(async()=>{const n=await import("./joraw2-Bb3_vNP4.js");if(typeof n.default!="function")throw new Error("JoRaw2 WASM import failed");const i=new URL("/assets/joraw2-3YkywkGx.wasm",import.meta.url).href;return n.default({locateFile(s,a){return s.endsWith("joraw2.wasm")?i:a+s}})})());const e=(await ot).LibRaw;if(!e)throw new Error("JoRaw2 class not found");if(!Rt){const n=new e;try{if(typeof n.runtimeInfo!="function")throw new Error("JoRaw2 runtime identity is missing");const i=n.runtimeInfo();if((i==null?void 0:i.wrapper)!=="joraw2"||!String((i==null?void 0:i.librawVersion)||"").startsWith("0.22.2")||!(i!=null&&i.nikonHe)||!(i!=null&&i.nikonHeStar))throw new Error(`Unexpected JoRaw2 runtime: ${JSON.stringify(i)}`);Rt=!0}finally{(t=n.delete)==null||t.call(n)}}return e}var Pn=`#!/usr/bin/env python3
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
`,An=`#!/usr/bin/env python3
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
`,In=`#!/usr/bin/env python3
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
`;const Lt=512,Tn=Cn["Sony ILCE-7RM5"].split(/\s+/).map(Number).filter(Number.isFinite);let lt=null;function Rn(r){if(r.byteLength<8)return null;const e=r.getUint16(0,!1);return e===18761?!0:e===19789?!1:null}function nn(r,e,t){const n=Math.min(r.length,e+t);let i="";for(let s=e;s<n;s++){const a=r[s];if(a===0)break;i+=String.fromCharCode(a)}return i.trim()}function K(r,e,t,n,i){const s=t===1||t===2||t===7?1:t===3||t===8?2:t===4||t===9?4:0;if(!s)return[];const a=s*n,o=a<=4?i:r.getUint32(i,e);if(o<0||o+a>r.byteLength)return[];const l=[];for(let h=0;h<n;h++){const c=o+h*s;t===1||t===2||t===7?l.push(r.getUint8(c)):t===3?l.push(r.getUint16(c,e)):t===8?l.push(r.getInt16(c,e)):t===4?l.push(r.getUint32(c,e)):t===9&&l.push(r.getInt32(c,e))}return l}function Ot(r,e,t,n,i,s){if(n!==2||i<=0)return"";const a=i<=4?s:e.getUint32(s,t);return a<0||a>=r.length?"":nn(r,a,i)}function rn(r){const e=new Uint8Array(r),t=new DataView(r),n=Rn(t);if(n===null||t.getUint16(2,n)!==42)return null;const s=d=>t.getUint16(d,n),a=d=>t.getUint32(d,n),o=[a(4)],l=new Set;let h="",c="";for(;o.length;){const d=o.pop();if(l.has(d)||d<=0||d+2>t.byteLength)continue;l.add(d);const u=s(d);if(d+2+u*12+4>t.byteLength)continue;const p=new Map;for(let T=0;T<u;T++){const P=d+2+T*12,O=s(P),E=s(P+2),U=a(P+4),k=P+8;p.set(O,{type:E,count:U,valueOffset:k})}const _=p.get(271),g=p.get(272);_&&!h&&(h=Ot(e,t,n,_.type,_.count,_.valueOffset)),g&&!c&&(c=Ot(e,t,n,g.type,g.count,g.valueOffset));const y=p.get(330);if(y){const T=K(t,n,y.type,y.count,y.valueOffset);for(const P of T)o.push(P)}const b=p.get(259),v=p.get(262);if(b&&v){const T=K(t,n,b.type,b.count,b.valueOffset)[0],P=K(t,n,v.type,v.count,v.valueOffset)[0];if(T===32766&&P===32803){const O=K(t,n,p.get(256).type,p.get(256).count,p.get(256).valueOffset)[0],E=K(t,n,p.get(257).type,p.get(257).count,p.get(257).valueOffset)[0],U=K(t,n,p.get(258).type,p.get(258).count,p.get(258).valueOffset)[0],k=K(t,n,p.get(273).type,p.get(273).count,p.get(273).valueOffset)[0],W=K(t,n,p.get(279).type,p.get(279).count,p.get(279).valueOffset)[0],q=p.get(33422)?K(t,n,p.get(33422).type,p.get(33422).count,p.get(33422).valueOffset):[0,1,1,2];p.get(29456)&&K(t,n,p.get(29456).type,p.get(29456).count,p.get(29456).valueOffset);const X=p.get(50717)?K(t,n,p.get(50717).type,p.get(50717).count,p.get(50717).valueOffset)[0]:16383,x=p.get(50719)?K(t,n,p.get(50719).type,p.get(50719).count,p.get(50719).valueOffset):[],B=p.get(50720)?K(t,n,p.get(50720).type,p.get(50720).count,p.get(50720).valueOffset):[];if(k+Lt+16>e.length||k+W>e.length)return null;const f=k+Lt,C=nn(e,f,4),D=e[f+8]<<8|e[f+9],S=e[f+10]<<8|e[f+11],w=e[f+12]<<8|e[f+13],R=e[f+14]<<8|e[f+15],J=w>>4&63,ue=R>>13,ne=R>>10&3,Q=S*2,V=W>=4?(e[k]|e[k+1]<<8|e[k+2]<<16|e[k+3]<<24)>>>0:0,nt=D===O&&Q===E;let ve=!1;if(V>=1&&V<=16&&W>=8+V*24){const re=new Map,se=new Map;let Te=!0;for(let Re=0;Re<V;Re++){const ae=k+8+Re*24,fe=t.getUint32(ae+8,!0),ge=t.getUint32(ae+12,!0),st=t.getUint32(ae+16,!0),Ee=t.getUint32(ae+20,!0);if(!st||!Ee||fe+st>O||ge+Ee>E){Te=!1;break}const Tt=re.get(ge);if(Tt!==void 0&&Tt!==Ee){Te=!1;break}re.set(ge,Ee),se.set(ge,(se.get(ge)||0)+st)}if(Te){const Re=Array.from(re.keys()).sort((fe,ge)=>fe-ge);let ae=0;for(const fe of Re){if(fe!==ae||se.get(fe)!==O){Te=!1;break}ae+=re.get(fe)}ve=Te&&ae===E}}const Me=V>=1&&V<=16&&D>0&&Q>0&&O%D===0&&E%Q===0&&V===O/D*(E/Q);if(C!=="A000"&&C!=="0000"||!nt&&!Me&&!ve||J!==16||ue!==3||ne!==3)return null;const rt=[1024,1024,1024,1024],it=q.length>=4?q.slice(0,4).map(re=>re===0?"R":re===2?"B":re===1?"G":"?").join(""):"";return{width:O,height:E,bitsPerSample:U,compression:T,photometric:P,blackLevel:rt,whiteLevel:Number(X||16383),cfaPattern:it,defaultCropOrigin:x.length>=2?[Number(x[0]),Number(x[1])]:void 0,defaultCropSize:B.length>=2?[Number(B[0]),Number(B[1])]:void 0,make:h||"SONY",model:c||"ILCE-7M5"}}}const A=a(d+2+u*12);A&&o.push(A)}return null}async function Ln(r){return lt||(lt=(async()=>{const e=await vn();return e.__jtrSonyCrawHqDecoderReady||(await e.FS.mkdirTree("/sony_craw_hq"),await e.FS.writeFile("/sony_craw_hq/llvc3_bitstream_probe.py",Pn),await e.FS.writeFile("/sony_craw_hq/llvc3_entropy.py",An),await e.FS.writeFile("/sony_craw_hq/llvc3_math.py",In),await e.runPythonAsync(`
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
`),e.__jtrSonyCrawHqDecoderReady=!0),e})()),lt}function On(r){return rn(r)}function Dn(r){return r==="ILCE-7M5"?Tn:null}function sn(r,e,t,n){var l;if(r.length!==e.width*e.height)throw new Error(`Sony cRAW HQ decoded size mismatch: got ${r.length}, expected ${e.width*e.height}`);const i=e.model||"ILCE-7M5",s=i.startsWith("Sony ")?i:`Sony ${i}`,a=Dn(i),o={...n||{},make:e.make||(n==null?void 0:n.camera_make)||"SONY",model:i,camera_make:e.make||(n==null?void 0:n.camera_make)||"SONY",camera_model:i,UniqueCameraModel:s,sourceFormat:t==="joraw2-wasm"?"Sony cRAW HQ / LLVC3 (JoRaw2 WASM)":"Sony cRAW HQ / LLVC3",sonyCrawHq:{...e,decodeBackend:t},color_desc:e.cfaPattern,black_level_per_channel:e.blackLevel,white_level:e.whiteLevel,color_matrix:a&&a.length===9?a:void 0,idata:{filters:2492765332,colors:3},color_data:{...(n==null?void 0:n.color_data)||{},black:1024,cblack_rawpy_style:e.blackLevel,dng_levels:{...((l=n==null?void 0:n.color_data)==null?void 0:l.dng_levels)||{},dng_cblack:e.blackLevel,dng_whitelevel:e.whiteLevel}}};return{data:r,width:e.width,height:e.height,bayerPattern:e.cfaPattern,blackLevels:e.blackLevel,whiteLevel:e.whiteLevel,metadata:o,isThreePlane:!1,isXTrans:!1}}async function Fn(r,e,t){const n=typeof performance<"u"?performance.now():Date.now(),i=await bt(),s=typeof performance<"u"?performance.now():Date.now(),a=new i;try{const o=new Uint8Array(r);await a.open(o,{});const l=typeof performance<"u"?performance.now():Date.now();let h=null;try{h=await a.metadata(!0)}catch(y){console.warn("[Sony cRAW HQ] fast WASM metadata read failed",y)}const c=typeof performance<"u"?performance.now():Date.now(),d=a.getRawImage(),u=typeof performance<"u"?performance.now():Date.now();if(!d||!d.data)throw new Error("Sony cRAW HQ LibRaw WASM returned no raw image");const p=d.data instanceof Uint16Array?d.data:new Uint16Array(d.data.buffer,d.data.byteOffset||0,d.data.byteLength/2),_=typeof performance<"u"?performance.now():Date.now();if(d.width!==e.width||d.height!==e.height)throw new Error(`Sony cRAW HQ LibRaw WASM dimensions mismatch: got ${d.width}x${d.height}, expected ${e.width}x${e.height}`);const g=typeof performance<"u"?performance.now():Date.now();return console.info("[Sony cRAW HQ] fast decode timings",{width:e.width,height:e.height,backend:"joraw2-wasm",wasmReadyMs:Math.round(s-n),openMs:Math.round(l-s),metadataMs:Math.round(c-l),unpackMs:Math.round(u-c),copyMs:Math.round(_-u),totalMs:Math.round(g-n)}),{rawImageData:sn(p,e,"joraw2-wasm",h),info:e}}finally{typeof a.delete=="function"?a.delete():typeof a.close=="function"&&a.close()}}async function Mn(r,e,t){const n=typeof performance<"u"?performance.now():Date.now(),i=await Ln(),s=typeof performance<"u"?performance.now():Date.now(),a=new Uint8Array(r),o=await fetch(new URL("/assets/sony_llvc3_static_lut4096_padded_u16-FsVBk-IV.bin",import.meta.url));if(!o.ok)throw new Error(`Failed to load Sony LLVC3 sample LUT: HTTP ${o.status}`);const l=new Uint8Array(await o.arrayBuffer()),h=typeof performance<"u"?performance.now():Date.now();i.globals.set("jtr_sony_arw_bytes",a),i.globals.set("jtr_sony_lut_bytes",l);const c=await i.runPythonAsync("jtr_decode_sony_craw_hq(jtr_sony_arw_bytes.to_py(), jtr_sony_lut_bytes.to_py())"),d=typeof performance<"u"?performance.now():Date.now(),u=c.toJs();typeof c.destroy=="function"&&c.destroy(),i.globals.delete("jtr_sony_arw_bytes"),i.globals.delete("jtr_sony_lut_bytes");const p=new Uint8Array(u.byteLength);p.set(u);const _=new Uint16Array(p.buffer),g=typeof performance<"u"?performance.now():Date.now(),y=typeof performance<"u"?performance.now():Date.now();return console.info("[Sony cRAW HQ] decode timings",{width:e.width,height:e.height,backend:"pyodide",pyodideReadyMs:Math.round(s-n),lutLoadMs:Math.round(h-s),llvc3DecodeMs:Math.round(d-h),copyMs:Math.round(g-d),totalMs:Math.round(y-n)}),{rawImageData:sn(_,e,"pyodide"),info:e}}async function En(r,e){const t=rn(r);if(!t)return null;try{return await Fn(r,t,e)}catch(n){return console.warn("[Sony cRAW HQ] fast WASM decode failed; falling back to Pyodide",n),Mn(r,t)}}async function Un(r,e){return En(r,e)}const Nn=["RGGB","BGGR","GRBG","GBRG"],Bn=new Set(Nn);function Vn(r){if(!r||typeof r!="object"||Array.isArray(r)||ArrayBuffer.isView(r))return r;const e=r;return e.value??e.values??e.description??r}function Ue(r){const e=Vn(r);let t="";if(typeof e=="string")t=e.toUpperCase().replace(/[^RGB012]/g,"");else if(typeof e=="number")t=String(e);else if(Array.isArray(e)||ArrayBuffer.isView(e))t=Array.from(e).map(String).join("");else return null;return/^[012]{4}$/.test(t)&&(t=t.replace(/0/g,"R").replace(/1/g,"G").replace(/2/g,"B")),Bn.has(t)?t:null}function Dt(r){var e;return r?[r.cfa_pattern,r.cfaPattern,r.BayerPattern,r.CFAPattern2,r.CFAPattern,(e=r.idata)==null?void 0:e.cfa_pattern]:[]}function zn(r,e,t){const n=Ue(r.bayerPattern),i=r.bayerPatternSource==="manual";if(n&&!i)return{pattern:n,source:r.bayerPatternSource||"decoder"};if(!i)for(const a of Dt(r.metadata)){const o=Ue(a);if(o)return{pattern:o,source:"libraw-metadata"}}for(const a of Dt(e)){const o=Ue(a);if(o)return{pattern:o,source:"metadata"}}const s=Ue(t)||(i?n:null);return s?{pattern:s,source:"manual"}:{pattern:null,source:null}}function Gn(r,e){return r.bayerPattern=e.pattern||"",r.bayerPatternSource=e.source||void 0,e.pattern&&e.source!=="manual"&&(r.metadata={...r.metadata||{},cfa_pattern:e.pattern,cfaPattern:e.pattern}),e}const Z=r=>{const e=Number(r);return Number.isFinite(e)?Math.max(0,e):0};function me(r){if(!r||typeof r.length!="number")return null;const e=Array.from(r);return e.length<4?null:[Z(e[0]),Z(e[1]),Z(e[2]),Z(e[3])]}function Ft(r){if(!r||r.source!=="libraw")return null;const e=me(r.channelOffsets),t=me(r.channelLevels),n=me(r.siteColorIndices),i=me(r.siteBaseLevels),s=me(r.siteLevels);if(!e||!t||!n||!i||!s)return null;const a=Math.max(0,Math.floor(Z(r.repeatRows))),o=Math.max(0,Math.floor(Z(r.repeatCols))),l=a*o,h=r.repeatValues&&typeof r.repeatValues.length=="number"?Array.from(r.repeatValues,Z):[],c=l>0&&l<=4098&&h.length>=l;return{source:"libraw",common:Z(r.common),channelOffsets:e,channelLevels:t,siteColorIndices:n,siteBaseLevels:i,siteLevels:s,repeatRows:c?a:0,repeatCols:c?o:0,repeatOriginY:Math.max(0,Math.floor(Z(r.repeatOriginY))),repeatOriginX:Math.max(0,Math.floor(Z(r.repeatOriginX))),repeatValues:c?h.slice(0,l):[]}}var $e=typeof self<"u"?self:global;const Fe=typeof navigator<"u",Hn=Fe&&typeof HTMLImageElement>"u",je=!(typeof global>"u"||typeof process>"u"||!process.versions||!process.versions.node),Je=$e.Buffer,Ne=$e.BigInt,Qe=!!Je,jn=r=>r;function We(r,e=jn){if(je)try{return typeof require=="function"?Promise.resolve(e(require(r))):import(r).then(e)}catch{console.warn(`Couldn't load ${r}`)}}let xt=$e.fetch;const Wn=r=>xt=r;if(!$e.fetch){const r=We("http",n=>n),e=We("https",n=>n),t=(n,{headers:i}={})=>new Promise(async(s,a)=>{let{port:o,hostname:l,pathname:h,protocol:c,search:d}=new URL(n);const u={method:"GET",hostname:l,path:encodeURI(h)+d,headers:i};o!==""&&(u.port=Number(o));const p=(c==="https:"?await e:await r).request(u,_=>{if(_.statusCode===301||_.statusCode===302){let g=new URL(_.headers.location,n).toString();return t(g,{headers:i}).then(s).catch(a)}s({status:_.statusCode,arrayBuffer:()=>new Promise(g=>{let y=[];_.on("data",b=>y.push(b)),_.on("end",()=>g(Buffer.concat(y)))})})});p.on("error",a),p.end()});Wn(t)}function m(r,e,t){return e in r?Object.defineProperty(r,e,{value:t,enumerable:!0,configurable:!0,writable:!0}):r[e]=t,r}const Xe=r=>an(r)?void 0:r,Xn=r=>r!==void 0;function an(r){return r===void 0||(r instanceof Map?r.size===0:Object.values(r).filter(Xn).length===0)}function M(r){let e=new Error(r);throw delete e.stack,e}function _e(r){return(r=function(e){for(;e.endsWith("\0");)e=e.slice(0,-1);return e}(r).trim())===""?void 0:r}function gt(r){let e=function(t){let n=0;return t.ifd0.enabled&&(n+=1024),t.exif.enabled&&(n+=2048),t.makerNote&&(n+=2048),t.userComment&&(n+=1024),t.gps.enabled&&(n+=512),t.interop.enabled&&(n+=100),t.ifd1.enabled&&(n+=1024),n+2048}(r);return r.jfif.enabled&&(e+=50),r.xmp.enabled&&(e+=2e4),r.iptc.enabled&&(e+=14e3),r.icc.enabled&&(e+=6e3),e}const mt=r=>String.fromCharCode.apply(null,r),Mt=typeof TextDecoder<"u"?new TextDecoder("utf-8"):void 0;function on(r){return Mt?Mt.decode(r):Qe?Buffer.from(r).toString("utf8"):decodeURIComponent(escape(mt(r)))}class j{static from(e,t){return e instanceof this&&e.le===t?e:new j(e,void 0,void 0,t)}constructor(e,t=0,n,i){if(typeof i=="boolean"&&(this.le=i),Array.isArray(e)&&(e=new Uint8Array(e)),e===0)this.byteOffset=0,this.byteLength=0;else if(e instanceof ArrayBuffer){n===void 0&&(n=e.byteLength-t);let s=new DataView(e,t,n);this._swapDataView(s)}else if(e instanceof Uint8Array||e instanceof DataView||e instanceof j){n===void 0&&(n=e.byteLength-t),(t+=e.byteOffset)+n>e.byteOffset+e.byteLength&&M("Creating view outside of available memory in ArrayBuffer");let s=new DataView(e.buffer,t,n);this._swapDataView(s)}else if(typeof e=="number"){let s=new DataView(new ArrayBuffer(e));this._swapDataView(s)}else M("Invalid input argument for BufferView: "+e)}_swapArrayBuffer(e){this._swapDataView(new DataView(e))}_swapBuffer(e){this._swapDataView(new DataView(e.buffer,e.byteOffset,e.byteLength))}_swapDataView(e){this.dataView=e,this.buffer=e.buffer,this.byteOffset=e.byteOffset,this.byteLength=e.byteLength}_lengthToEnd(e){return this.byteLength-e}set(e,t,n=j){return e instanceof DataView||e instanceof j?e=new Uint8Array(e.buffer,e.byteOffset,e.byteLength):e instanceof ArrayBuffer&&(e=new Uint8Array(e)),e instanceof Uint8Array||M("BufferView.set(): Invalid data argument."),this.toUint8().set(e,t),new n(this,t,e.byteLength)}subarray(e,t){return t=t||this._lengthToEnd(e),new j(this,e,t)}toUint8(){return new Uint8Array(this.buffer,this.byteOffset,this.byteLength)}getUint8Array(e,t){return new Uint8Array(this.buffer,this.byteOffset+e,t)}getString(e=0,t=this.byteLength){return on(this.getUint8Array(e,t))}getLatin1String(e=0,t=this.byteLength){let n=this.getUint8Array(e,t);return mt(n)}getUnicodeString(e=0,t=this.byteLength){const n=[];for(let i=0;i<t&&e+i<this.byteLength;i+=2)n.push(this.getUint16(e+i));return mt(n)}getInt8(e){return this.dataView.getInt8(e)}getUint8(e){return this.dataView.getUint8(e)}getInt16(e,t=this.le){return this.dataView.getInt16(e,t)}getInt32(e,t=this.le){return this.dataView.getInt32(e,t)}getUint16(e,t=this.le){return this.dataView.getUint16(e,t)}getUint32(e,t=this.le){return this.dataView.getUint32(e,t)}getFloat32(e,t=this.le){return this.dataView.getFloat32(e,t)}getFloat64(e,t=this.le){return this.dataView.getFloat64(e,t)}getFloat(e,t=this.le){return this.dataView.getFloat32(e,t)}getDouble(e,t=this.le){return this.dataView.getFloat64(e,t)}getUintBytes(e,t,n){switch(t){case 1:return this.getUint8(e,n);case 2:return this.getUint16(e,n);case 4:return this.getUint32(e,n);case 8:return this.getUint64&&this.getUint64(e,n)}}getUint(e,t,n){switch(t){case 8:return this.getUint8(e,n);case 16:return this.getUint16(e,n);case 32:return this.getUint32(e,n);case 64:return this.getUint64&&this.getUint64(e,n)}}toString(e){return this.dataView.toString(e,this.constructor.name)}ensureChunk(){}}function _t(r,e){M(`${r} '${e}' was not loaded, try using full build of exifr.`)}class St extends Map{constructor(e){super(),this.kind=e}get(e,t){return this.has(e)||_t(this.kind,e),t&&(e in t||function(n,i){M(`Unknown ${n} '${i}'.`)}(this.kind,e),t[e].enabled||_t(this.kind,e)),super.get(e)}keyList(){return Array.from(this.keys())}}var ee=new St("file parser"),F=new St("segment parser"),te=new St("file reader");function Kn(r,e){return typeof r=="string"?Et(r,e):Fe&&!Hn&&r instanceof HTMLImageElement?Et(r.src,e):r instanceof Uint8Array||r instanceof ArrayBuffer||r instanceof DataView?new j(r):Fe&&r instanceof Blob?yt(r,e,"blob",xe):void M("Invalid input argument")}function Et(r,e){return(t=r).startsWith("data:")||t.length>1e4?wt(r,e,"base64"):je&&r.includes("://")?yt(r,e,"url",be):je?wt(r,e,"fs"):Fe?yt(r,e,"url",be):void M("Invalid input argument");var t}async function yt(r,e,t,n){return te.has(t)?wt(r,e,t):n?async function(i,s){let a=await s(i);return new j(a)}(r,n):void M(`Parser ${t} is not loaded`)}async function wt(r,e,t){let n=new(te.get(t))(r,e);return await n.read(),n}const be=r=>xt(r).then(e=>e.arrayBuffer()),xe=r=>new Promise((e,t)=>{let n=new FileReader;n.onloadend=()=>e(n.result||new ArrayBuffer),n.onerror=t,n.readAsArrayBuffer(r)});class qn extends Map{get tagKeys(){return this.allKeys||(this.allKeys=Array.from(this.keys())),this.allKeys}get tagValues(){return this.allValues||(this.allValues=Array.from(this.values())),this.allValues}}function L(r,e,t){let n=new qn;for(let[i,s]of t)n.set(i,s);if(Array.isArray(e))for(let i of e)r.set(i,n);else r.set(e,n);return n}function Se(r,e,t){let n,i=r.get(e);for(n of t)i.set(n[0],n[1])}const N=new Map,$=new Map,he=new Map,oe=["chunked","firstChunkSize","firstChunkSizeNode","firstChunkSizeBrowser","chunkSize","chunkLimit"],Ae=["jfif","xmp","icc","iptc","ihdr"],ke=["tiff",...Ae],I=["ifd0","ifd1","exif","gps","interop"],le=[...ke,...I],de=["makerNote","userComment"],Ie=["translateKeys","translateValues","reviveValues","multiSegment"],ce=[...Ie,"sanitize","mergeOutput","silentErrors"];class ln{get translate(){return this.translateKeys||this.translateValues||this.reviveValues}}class Le extends ln{get needed(){return this.enabled||this.deps.size>0}constructor(e,t,n,i){if(super(),m(this,"enabled",!1),m(this,"skip",new Set),m(this,"pick",new Set),m(this,"deps",new Set),m(this,"translateKeys",!1),m(this,"translateValues",!1),m(this,"reviveValues",!1),this.key=e,this.enabled=t,this.parse=this.enabled,this.applyInheritables(i),this.canBeFiltered=I.includes(e),this.canBeFiltered&&(this.dict=N.get(e)),n!==void 0)if(Array.isArray(n))this.parse=this.enabled=!0,this.canBeFiltered&&n.length>0&&this.translateTagSet(n,this.pick);else if(typeof n=="object"){if(this.enabled=!0,this.parse=n.parse!==!1,this.canBeFiltered){let{pick:s,skip:a}=n;s&&s.length>0&&this.translateTagSet(s,this.pick),a&&a.length>0&&this.translateTagSet(a,this.skip)}this.applyInheritables(n)}else n===!0||n===!1?this.parse=this.enabled=n:M(`Invalid options argument: ${n}`)}applyInheritables(e){let t,n;for(t of Ie)n=e[t],n!==void 0&&(this[t]=n)}translateTagSet(e,t){if(this.dict){let n,i,{tagKeys:s,tagValues:a}=this.dict;for(n of e)typeof n=="string"?(i=a.indexOf(n),i===-1&&(i=s.indexOf(Number(n))),i!==-1&&t.add(Number(s[i]))):t.add(n)}else for(let n of e)t.add(n)}finalizeFilters(){!this.enabled&&this.deps.size>0?(this.enabled=!0,Ke(this.pick,this.deps)):this.enabled&&this.pick.size>0&&Ke(this.pick,this.deps)}}var G={jfif:!1,tiff:!0,xmp:!1,icc:!1,iptc:!1,ifd0:!0,ifd1:!1,exif:!0,gps:!0,interop:!1,ihdr:void 0,makerNote:!1,userComment:!1,multiSegment:!1,skip:[],pick:[],translateKeys:!0,translateValues:!0,reviveValues:!0,sanitize:!0,mergeOutput:!0,silentErrors:!0,chunked:!0,firstChunkSize:void 0,firstChunkSizeNode:512,firstChunkSizeBrowser:65536,chunkSize:65536,chunkLimit:5},Ut=new Map;class Ce extends ln{static useCached(e){let t=Ut.get(e);return t!==void 0||(t=new this(e),Ut.set(e,t)),t}constructor(e){super(),e===!0?this.setupFromTrue():e===void 0?this.setupFromUndefined():Array.isArray(e)?this.setupFromArray(e):typeof e=="object"?this.setupFromObject(e):M(`Invalid options argument ${e}`),this.firstChunkSize===void 0&&(this.firstChunkSize=Fe?this.firstChunkSizeBrowser:this.firstChunkSizeNode),this.mergeOutput&&(this.ifd1.enabled=!1),this.filterNestedSegmentTags(),this.traverseTiffDependencyTree(),this.checkLoadedPlugins()}setupFromUndefined(){let e;for(e of oe)this[e]=G[e];for(e of ce)this[e]=G[e];for(e of de)this[e]=G[e];for(e of le)this[e]=new Le(e,G[e],void 0,this)}setupFromTrue(){let e;for(e of oe)this[e]=G[e];for(e of ce)this[e]=G[e];for(e of de)this[e]=!0;for(e of le)this[e]=new Le(e,!0,void 0,this)}setupFromArray(e){let t;for(t of oe)this[t]=G[t];for(t of ce)this[t]=G[t];for(t of de)this[t]=G[t];for(t of le)this[t]=new Le(t,!1,void 0,this);this.setupGlobalFilters(e,void 0,I)}setupFromObject(e){let t;for(t of(I.ifd0=I.ifd0||I.image,I.ifd1=I.ifd1||I.thumbnail,Object.assign(this,e),oe))this[t]=dt(e[t],G[t]);for(t of ce)this[t]=dt(e[t],G[t]);for(t of de)this[t]=dt(e[t],G[t]);for(t of ke)this[t]=new Le(t,G[t],e[t],this);for(t of I)this[t]=new Le(t,G[t],e[t],this.tiff);this.setupGlobalFilters(e.pick,e.skip,I,le),e.tiff===!0?this.batchEnableWithBool(I,!0):e.tiff===!1?this.batchEnableWithUserValue(I,e):Array.isArray(e.tiff)?this.setupGlobalFilters(e.tiff,void 0,I):typeof e.tiff=="object"&&this.setupGlobalFilters(e.tiff.pick,e.tiff.skip,I)}batchEnableWithBool(e,t){for(let n of e)this[n].enabled=t}batchEnableWithUserValue(e,t){for(let n of e){let i=t[n];this[n].enabled=i!==!1&&i!==void 0}}setupGlobalFilters(e,t,n,i=n){if(e&&e.length){for(let a of i)this[a].enabled=!1;let s=Nt(e,n);for(let[a,o]of s)Ke(this[a].pick,o),this[a].enabled=!0}else if(t&&t.length){let s=Nt(t,n);for(let[a,o]of s)Ke(this[a].skip,o)}}filterNestedSegmentTags(){let{ifd0:e,exif:t,xmp:n,iptc:i,icc:s}=this;this.makerNote?t.deps.add(37500):t.skip.add(37500),this.userComment?t.deps.add(37510):t.skip.add(37510),n.enabled||e.skip.add(700),i.enabled||e.skip.add(33723),s.enabled||e.skip.add(34675)}traverseTiffDependencyTree(){let{ifd0:e,exif:t,gps:n,interop:i}=this;i.needed&&(t.deps.add(40965),e.deps.add(40965)),t.needed&&e.deps.add(34665),n.needed&&e.deps.add(34853),this.tiff.enabled=I.some(s=>this[s].enabled===!0)||this.makerNote||this.userComment;for(let s of I)this[s].finalizeFilters()}get onlyTiff(){return!Ae.map(e=>this[e].enabled).some(e=>e===!0)&&this.tiff.enabled}checkLoadedPlugins(){for(let e of ke)this[e].enabled&&!F.has(e)&&_t("segment parser",e)}}function Nt(r,e){let t,n,i,s,a=[];for(i of e){for(s of(t=N.get(i),n=[],t))(r.includes(s[0])||r.includes(s[1]))&&n.push(s[0]);n.length&&a.push([i,n])}return a}function dt(r,e){return r!==void 0?r:e!==void 0?e:void 0}function Ke(r,e){for(let t of e)r.add(t)}m(Ce,"default",G);class pe{constructor(e){m(this,"parsers",{}),m(this,"output",{}),m(this,"errors",[]),m(this,"pushToErrors",t=>this.errors.push(t)),this.options=Ce.useCached(e)}async read(e){this.file=await Kn(e,this.options)}setup(){if(this.fileParser)return;let{file:e}=this,t=e.getUint16(0);for(let[n,i]of ee)if(i.canHandle(e,t))return this.fileParser=new i(this.options,this.file,this.parsers),e[n]=!0;this.file.close&&this.file.close(),M("Unknown file format")}async parse(){let{output:e,errors:t}=this;return this.setup(),this.options.silentErrors?(await this.executeParsers().catch(this.pushToErrors),t.push(...this.fileParser.errors)):await this.executeParsers(),this.file.close&&this.file.close(),this.options.silentErrors&&t.length>0&&(e.errors=t),Xe(e)}async executeParsers(){let{output:e}=this;await this.fileParser.parse();let t=Object.values(this.parsers).map(async n=>{let i=await n.parse();n.assignToOutput(e,i)});this.options.silentErrors&&(t=t.map(n=>n.catch(this.pushToErrors))),await Promise.all(t)}async extractThumbnail(){this.setup();let{options:e,file:t}=this,n=F.get("tiff",e);var i;if(t.tiff?i={start:0,type:"tiff"}:t.jpeg&&(i=await this.fileParser.getOrFindSegment("tiff")),i===void 0)return;let s=await this.fileParser.ensureSegmentChunk(i),a=this.parsers.tiff=new n(s,e,t),o=await a.extractThumbnail();return t.close&&t.close(),o}}async function Ze(r,e){let t=new pe(e);return await t.read(r),t.parse()}var Yn=Object.freeze({__proto__:null,parse:Ze,Exifr:pe,fileParsers:ee,segmentParsers:F,fileReaders:te,tagKeys:N,tagValues:$,tagRevivers:he,createDictionary:L,extendDictionary:Se,fetchUrlAsArrayBuffer:be,readBlobAsArrayBuffer:xe,chunkedProps:oe,otherSegments:Ae,segments:ke,tiffBlocks:I,segmentsAndBlocks:le,tiffExtractables:de,inheritables:Ie,allFormatters:ce,Options:Ce});class et{constructor(e,t,n){m(this,"errors",[]),m(this,"ensureSegmentChunk",async i=>{let s=i.start,a=i.size||65536;if(this.file.chunked)if(this.file.available(s,a))i.chunk=this.file.subarray(s,a);else try{i.chunk=await this.file.readChunk(s,a)}catch(o){M(`Couldn't read segment: ${JSON.stringify(i)}. ${o.message}`)}else this.file.byteLength>s+a?i.chunk=this.file.subarray(s,a):i.size===void 0?i.chunk=this.file.subarray(s):M("Segment unreachable: "+JSON.stringify(i));return i.chunk}),this.extendOptions&&this.extendOptions(e),this.options=e,this.file=t,this.parsers=n}injectSegment(e,t){this.options[e].enabled&&this.createParser(e,t)}createParser(e,t){let n=new(F.get(e))(t,this.options,this.file);return this.parsers[e]=n}createParsers(e){for(let t of e){let{type:n,chunk:i}=t,s=this.options[n];if(s&&s.enabled){let a=this.parsers[n];a&&a.append||a||this.createParser(n,i)}}}async readSegments(e){let t=e.map(this.ensureSegmentChunk);await Promise.all(t)}}class Y{static findPosition(e,t){let n=e.getUint16(t+2)+2,i=typeof this.headerLength=="function"?this.headerLength(e,t,n):this.headerLength,s=t+i,a=n-i;return{offset:t,length:n,headerLength:i,start:s,size:a,end:s+a}}static parse(e,t={}){return new this(e,new Ce({[this.type]:t}),e).parse()}normalizeInput(e){return e instanceof j?e:new j(e)}constructor(e,t={},n){m(this,"errors",[]),m(this,"raw",new Map),m(this,"handleError",i=>{if(!this.options.silentErrors)throw i;this.errors.push(i.message)}),this.chunk=this.normalizeInput(e),this.file=n,this.type=this.constructor.type,this.globalOptions=this.options=t,this.localOptions=t[this.type],this.canTranslate=this.localOptions&&this.localOptions.translate}translate(){this.canTranslate&&(this.translated=this.translateBlock(this.raw,this.type))}get output(){return this.translated?this.translated:this.raw?Object.fromEntries(this.raw):void 0}translateBlock(e,t){let n=he.get(t),i=$.get(t),s=N.get(t),a=this.options[t],o=a.reviveValues&&!!n,l=a.translateValues&&!!i,h=a.translateKeys&&!!s,c={};for(let[d,u]of e)o&&n.has(d)?u=n.get(d)(u):l&&i.has(d)&&(u=this.translateValue(u,i.get(d))),h&&s.has(d)&&(d=s.get(d)||d),c[d]=u;return c}translateValue(e,t){return t[e]||t.DEFAULT||e}assignToOutput(e,t){this.assignObjectToOutput(e,this.constructor.type,t)}assignObjectToOutput(e,t,n){if(this.globalOptions.mergeOutput)return Object.assign(e,n);e[t]?Object.assign(e[t],n):e[t]=n}}m(Y,"headerLength",4),m(Y,"type",void 0),m(Y,"multiSegment",!1),m(Y,"canHandle",()=>!1);function $n(r){return r===192||r===194||r===196||r===219||r===221||r===218||r===254}function Jn(r){return r>=224&&r<=239}function Qn(r,e,t){for(let[n,i]of F)if(i.canHandle(r,e,t))return n}class Bt extends et{constructor(...e){super(...e),m(this,"appSegments",[]),m(this,"jpegSegments",[]),m(this,"unknownSegments",[])}static canHandle(e,t){return t===65496}async parse(){await this.findAppSegments(),await this.readSegments(this.appSegments),this.mergeMultiSegments(),this.createParsers(this.mergedAppSegments||this.appSegments)}setupSegmentFinderArgs(e){e===!0?(this.findAll=!0,this.wanted=new Set(F.keyList())):(e=e===void 0?F.keyList().filter(t=>this.options[t].enabled):e.filter(t=>this.options[t].enabled&&F.has(t)),this.findAll=!1,this.remaining=new Set(e),this.wanted=new Set(e)),this.unfinishedMultiSegment=!1}async findAppSegments(e=0,t){this.setupSegmentFinderArgs(t);let{file:n,findAll:i,wanted:s,remaining:a}=this;if(!i&&this.file.chunked&&(i=Array.from(s).some(o=>{let l=F.get(o),h=this.options[o];return l.multiSegment&&h.multiSegment}),i&&await this.file.readWhole()),e=this.findAppSegmentsInRange(e,n.byteLength),!this.options.onlyTiff&&n.chunked){let o=!1;for(;a.size>0&&!o&&(n.canReadNextChunk||this.unfinishedMultiSegment);){let{nextChunkOffset:l}=n,h=this.appSegments.some(c=>!this.file.available(c.offset||c.start,c.length||c.size));if(o=e>l&&!h?!await n.readNextChunk(e):!await n.readNextChunk(l),(e=this.findAppSegmentsInRange(e,n.byteLength))===void 0)return}}}findAppSegmentsInRange(e,t){t-=2;let n,i,s,a,o,l,{file:h,findAll:c,wanted:d,remaining:u,options:p}=this;for(;e<t;e++)if(h.getUint8(e)===255){if(n=h.getUint8(e+1),Jn(n)){if(i=h.getUint16(e+2),s=Qn(h,e,i),s&&d.has(s)&&(a=F.get(s),o=a.findPosition(h,e),l=p[s],o.type=s,this.appSegments.push(o),!c&&(a.multiSegment&&l.multiSegment?(this.unfinishedMultiSegment=o.chunkNumber<o.chunkCount,this.unfinishedMultiSegment||u.delete(s)):u.delete(s),u.size===0)))break;p.recordUnknownSegments&&(o=Y.findPosition(h,e),o.marker=n,this.unknownSegments.push(o)),e+=i+1}else if($n(n)){if(i=h.getUint16(e+2),n===218&&p.stopAfterSos!==!1)return;p.recordJpegSegments&&this.jpegSegments.push({offset:e,length:i,marker:n}),e+=i+1}}return e}mergeMultiSegments(){if(!this.appSegments.some(t=>t.multiSegment))return;let e=function(t,n){let i,s,a,o=new Map;for(let l=0;l<t.length;l++)i=t[l],s=i[n],o.has(s)?a=o.get(s):o.set(s,a=[]),a.push(i);return Array.from(o)}(this.appSegments,"type");this.mergedAppSegments=e.map(([t,n])=>{let i=F.get(t,this.options);return i.handleMultiSegments?{type:t,chunk:i.handleMultiSegments(n)}:n[0]})}getSegment(e){return this.appSegments.find(t=>t.type===e)}async getOrFindSegment(e){let t=this.getSegment(e);return t===void 0&&(await this.findAppSegments(0,[e]),t=this.getSegment(e)),t}}m(Bt,"type","jpeg"),ee.set("jpeg",Bt);const Zn=[void 0,1,1,2,4,8,1,1,2,4,8,4,8,4];class er extends Y{parseHeader(){var e=this.chunk.getUint16();e===18761?this.le=!0:e===19789&&(this.le=!1),this.chunk.le=this.le,this.headerParsed=!0}parseTags(e,t,n=new Map){let{pick:i,skip:s}=this.options[t];i=new Set(i);let a=i.size>0,o=s.size===0,l=this.chunk.getUint16(e);e+=2;for(let h=0;h<l;h++){let c=this.chunk.getUint16(e);if(a){if(i.has(c)&&(n.set(c,this.parseTag(e,c,t)),i.delete(c),i.size===0))break}else!o&&s.has(c)||n.set(c,this.parseTag(e,c,t));e+=12}return n}parseTag(e,t,n){let{chunk:i}=this,s=i.getUint16(e+2),a=i.getUint32(e+4),o=Zn[s];if(o*a<=4?e+=8:e=i.getUint32(e+8),(s<1||s>13)&&M(`Invalid TIFF value type. block: ${n.toUpperCase()}, tag: ${t.toString(16)}, type: ${s}, offset ${e}`),e>i.byteLength&&M(`Invalid TIFF value offset. block: ${n.toUpperCase()}, tag: ${t.toString(16)}, type: ${s}, offset ${e} is outside of chunk size ${i.byteLength}`),s===1)return i.getUint8Array(e,a);if(s===2)return _e(i.getString(e,a));if(s===7)return i.getUint8Array(e,a);if(a===1)return this.parseTagValue(s,e);{let l=new(function(c){switch(c){case 1:return Uint8Array;case 3:return Uint16Array;case 4:return Uint32Array;case 5:return Array;case 6:return Int8Array;case 8:return Int16Array;case 9:return Int32Array;case 10:return Array;case 11:return Float32Array;case 12:return Float64Array;default:return Array}}(s))(a),h=o;for(let c=0;c<a;c++)l[c]=this.parseTagValue(s,e),e+=h;return l}}parseTagValue(e,t){let{chunk:n}=this;switch(e){case 1:return n.getUint8(t);case 3:return n.getUint16(t);case 4:return n.getUint32(t);case 5:return n.getUint32(t)/n.getUint32(t+4);case 6:return n.getInt8(t);case 8:return n.getInt16(t);case 9:return n.getInt32(t);case 10:return n.getInt32(t)/n.getInt32(t+4);case 11:return n.getFloat(t);case 12:return n.getDouble(t);case 13:return n.getUint32(t);default:M(`Invalid tiff type ${e}`)}}}class ct extends er{static canHandle(e,t){return e.getUint8(t+1)===225&&e.getUint32(t+4)===1165519206&&e.getUint16(t+8)===0}async parse(){this.parseHeader();let{options:e}=this;return e.ifd0.enabled&&await this.parseIfd0Block(),e.exif.enabled&&await this.safeParse("parseExifBlock"),e.gps.enabled&&await this.safeParse("parseGpsBlock"),e.interop.enabled&&await this.safeParse("parseInteropBlock"),e.ifd1.enabled&&await this.safeParse("parseThumbnailBlock"),this.createOutput()}safeParse(e){let t=this[e]();return t.catch!==void 0&&(t=t.catch(this.handleError)),t}findIfd0Offset(){this.ifd0Offset===void 0&&(this.ifd0Offset=this.chunk.getUint32(4))}findIfd1Offset(){if(this.ifd1Offset===void 0){this.findIfd0Offset();let e=this.chunk.getUint16(this.ifd0Offset),t=this.ifd0Offset+2+12*e;this.ifd1Offset=this.chunk.getUint32(t)}}parseBlock(e,t){let n=new Map;return this[t]=n,this.parseTags(e,t,n),n}async parseIfd0Block(){if(this.ifd0)return;let{file:e}=this;this.findIfd0Offset(),this.ifd0Offset<8&&M("Malformed EXIF data"),!e.chunked&&this.ifd0Offset>e.byteLength&&M(`IFD0 offset points to outside of file.
this.ifd0Offset: ${this.ifd0Offset}, file.byteLength: ${e.byteLength}`),e.tiff&&await e.ensureChunk(this.ifd0Offset,gt(this.options));let t=this.parseBlock(this.ifd0Offset,"ifd0");return t.size!==0?(this.exifOffset=t.get(34665),this.interopOffset=t.get(40965),this.gpsOffset=t.get(34853),this.xmp=t.get(700),this.iptc=t.get(33723),this.icc=t.get(34675),this.options.sanitize&&(t.delete(34665),t.delete(40965),t.delete(34853),t.delete(700),t.delete(33723),t.delete(34675)),t):void 0}async parseExifBlock(){if(this.exif||(this.ifd0||await this.parseIfd0Block(),this.exifOffset===void 0))return;this.file.tiff&&await this.file.ensureChunk(this.exifOffset,gt(this.options));let e=this.parseBlock(this.exifOffset,"exif");return this.interopOffset||(this.interopOffset=e.get(40965)),this.makerNote=e.get(37500),this.userComment=e.get(37510),this.options.sanitize&&(e.delete(40965),e.delete(37500),e.delete(37510)),this.unpack(e,41728),this.unpack(e,41729),e}unpack(e,t){let n=e.get(t);n&&n.length===1&&e.set(t,n[0])}async parseGpsBlock(){if(this.gps||(this.ifd0||await this.parseIfd0Block(),this.gpsOffset===void 0))return;let e=this.parseBlock(this.gpsOffset,"gps");return e&&e.has(2)&&e.has(4)&&(e.set("latitude",Vt(...e.get(2),e.get(1))),e.set("longitude",Vt(...e.get(4),e.get(3)))),e}async parseInteropBlock(){if(!this.interop&&(this.ifd0||await this.parseIfd0Block(),this.interopOffset!==void 0||this.exif||await this.parseExifBlock(),this.interopOffset!==void 0))return this.parseBlock(this.interopOffset,"interop")}async parseThumbnailBlock(e=!1){if(!this.ifd1&&!this.ifd1Parsed&&(!this.options.mergeOutput||e))return this.findIfd1Offset(),this.ifd1Offset>0&&(this.parseBlock(this.ifd1Offset,"ifd1"),this.ifd1Parsed=!0),this.ifd1}async extractThumbnail(){if(this.headerParsed||this.parseHeader(),this.ifd1Parsed||await this.parseThumbnailBlock(!0),this.ifd1===void 0)return;let e=this.ifd1.get(513),t=this.ifd1.get(514);return this.chunk.getUint8Array(e,t)}get image(){return this.ifd0}get thumbnail(){return this.ifd1}createOutput(){let e,t,n,i={};for(t of I)if(e=this[t],!an(e))if(n=this.canTranslate?this.translateBlock(e,t):Object.fromEntries(e),this.options.mergeOutput){if(t==="ifd1")continue;Object.assign(i,n)}else i[t]=n;return this.makerNote&&(i.makerNote=this.makerNote),this.userComment&&(i.userComment=this.userComment),i}assignToOutput(e,t){if(this.globalOptions.mergeOutput)Object.assign(e,t);else for(let[n,i]of Object.entries(t))this.assignObjectToOutput(e,n,i)}}function Vt(r,e,t,n){var i=r+e/60+t/3600;return n!=="S"&&n!=="W"||(i*=-1),i}m(ct,"type","tiff"),m(ct,"headerLength",10),F.set("tiff",ct);var tr=Object.freeze({__proto__:null,default:Yn,Exifr:pe,fileParsers:ee,segmentParsers:F,fileReaders:te,tagKeys:N,tagValues:$,tagRevivers:he,createDictionary:L,extendDictionary:Se,fetchUrlAsArrayBuffer:be,readBlobAsArrayBuffer:xe,chunkedProps:oe,otherSegments:Ae,segments:ke,tiffBlocks:I,segmentsAndBlocks:le,tiffExtractables:de,inheritables:Ie,allFormatters:ce,Options:Ce,parse:Ze});const kt={ifd0:!1,ifd1:!1,exif:!1,gps:!1,interop:!1,sanitize:!1,reviveValues:!0,translateKeys:!1,translateValues:!1,mergeOutput:!1},Ct=Object.assign({},kt,{firstChunkSize:4e4,gps:[1,2,3,4]});async function dn(r){let e=new pe(Ct);await e.read(r);let t=await e.parse();if(t&&t.gps){let{latitude:n,longitude:i}=t.gps;return{latitude:n,longitude:i}}}const vt=Object.assign({},kt,{tiff:!1,ifd1:!0,mergeOutput:!1});async function cn(r){let e=new pe(vt);await e.read(r);let t=await e.extractThumbnail();return t&&Qe?Je.from(t):t}async function hn(r){let e=await this.thumbnail(r);if(e!==void 0){let t=new Blob([e]);return URL.createObjectURL(t)}}const Pt=Object.assign({},kt,{firstChunkSize:4e4,ifd0:[274]});async function At(r){let e=new pe(Pt);await e.read(r);let t=await e.parse();if(t&&t.ifd0)return t.ifd0[274]}const It=Object.freeze({1:{dimensionSwapped:!1,scaleX:1,scaleY:1,deg:0,rad:0},2:{dimensionSwapped:!1,scaleX:-1,scaleY:1,deg:0,rad:0},3:{dimensionSwapped:!1,scaleX:1,scaleY:1,deg:180,rad:180*Math.PI/180},4:{dimensionSwapped:!1,scaleX:-1,scaleY:1,deg:180,rad:180*Math.PI/180},5:{dimensionSwapped:!0,scaleX:1,scaleY:-1,deg:90,rad:90*Math.PI/180},6:{dimensionSwapped:!0,scaleX:1,scaleY:1,deg:90,rad:90*Math.PI/180},7:{dimensionSwapped:!0,scaleX:1,scaleY:-1,deg:270,rad:270*Math.PI/180},8:{dimensionSwapped:!0,scaleX:1,scaleY:1,deg:270,rad:270*Math.PI/180}});let ye=!0,we=!0;if(typeof navigator=="object"){let r=navigator.userAgent;if(r.includes("iPad")||r.includes("iPhone")){let e=r.match(/OS (\d+)_(\d+)/);if(e){let[,t,n]=e;ye=Number(t)+.1*Number(n)<13.4,we=!1}}else if(r.includes("OS X 10")){let[,e]=r.match(/OS X 10[_.](\d+)/);ye=we=Number(e)<15}if(r.includes("Chrome/")){let[,e]=r.match(/Chrome\/(\d+)/);ye=we=Number(e)<81}else if(r.includes("Firefox/")){let[,e]=r.match(/Firefox\/(\d+)/);ye=we=Number(e)<77}}async function pn(r){let e=await At(r);return Object.assign({canvas:ye,css:we},It[e])}class nr extends j{constructor(...e){super(...e),m(this,"ranges",new rr),this.byteLength!==0&&this.ranges.add(0,this.byteLength)}_tryExtend(e,t,n){if(e===0&&this.byteLength===0&&n){let i=new DataView(n.buffer||n,n.byteOffset,n.byteLength);this._swapDataView(i)}else{let i=e+t;if(i>this.byteLength){let{dataView:s}=this._extend(i);this._swapDataView(s)}}}_extend(e){let t;t=Qe?Je.allocUnsafe(e):new Uint8Array(e);let n=new DataView(t.buffer,t.byteOffset,t.byteLength);return t.set(new Uint8Array(this.buffer,this.byteOffset,this.byteLength),0),{uintView:t,dataView:n}}subarray(e,t,n=!1){return t=t||this._lengthToEnd(e),n&&this._tryExtend(e,t),this.ranges.add(e,t),super.subarray(e,t)}set(e,t,n=!1){n&&this._tryExtend(t,e.byteLength,e);let i=super.set(e,t);return this.ranges.add(t,i.byteLength),i}async ensureChunk(e,t){this.chunked&&(this.ranges.available(e,t)||await this.readChunk(e,t))}available(e,t){return this.ranges.available(e,t)}}class rr{constructor(){m(this,"list",[])}get length(){return this.list.length}add(e,t,n=0){let i=e+t,s=this.list.filter(a=>zt(e,a.offset,i)||zt(e,a.end,i));if(s.length>0){e=Math.min(e,...s.map(o=>o.offset)),i=Math.max(i,...s.map(o=>o.end)),t=i-e;let a=s.shift();a.offset=e,a.length=t,a.end=i,this.list=this.list.filter(o=>!s.includes(o))}else this.list.push({offset:e,length:t,end:i})}available(e,t){let n=e+t;return this.list.some(i=>i.offset<=e&&n<=i.end)}}function zt(r,e,t){return r<=e&&e<=t}class tt extends nr{constructor(e,t){super(0),m(this,"chunksRead",0),this.input=e,this.options=t}async readWhole(){this.chunked=!1,await this.readChunk(this.nextChunkOffset)}async readChunked(){this.chunked=!0,await this.readChunk(0,this.options.firstChunkSize)}async readNextChunk(e=this.nextChunkOffset){if(this.fullyRead)return this.chunksRead++,!1;let t=this.options.chunkSize,n=await this.readChunk(e,t);return!!n&&n.byteLength===t}async readChunk(e,t){if(this.chunksRead++,(t=this.safeWrapAddress(e,t))!==0)return this._readChunk(e,t)}safeWrapAddress(e,t){return this.size!==void 0&&e+t>this.size?Math.max(0,this.size-e):t}get nextChunkOffset(){if(this.ranges.list.length!==0)return this.ranges.list[0].length}get canReadNextChunk(){return this.chunksRead<this.options.chunkLimit}get fullyRead(){return this.size!==void 0&&this.nextChunkOffset===this.size}read(){return this.options.chunked?this.readChunked():this.readWhole()}close(){}}te.set("blob",class extends tt{async readWhole(){this.chunked=!1;let r=await xe(this.input);this._swapArrayBuffer(r)}readChunked(){return this.chunked=!0,this.size=this.input.size,super.readChunked()}async _readChunk(r,e){let t=e?r+e:void 0,n=this.input.slice(r,t),i=await xe(n);return this.set(i,r,!0)}});var ir=Object.freeze({__proto__:null,default:tr,Exifr:pe,fileParsers:ee,segmentParsers:F,fileReaders:te,tagKeys:N,tagValues:$,tagRevivers:he,createDictionary:L,extendDictionary:Se,fetchUrlAsArrayBuffer:be,readBlobAsArrayBuffer:xe,chunkedProps:oe,otherSegments:Ae,segments:ke,tiffBlocks:I,segmentsAndBlocks:le,tiffExtractables:de,inheritables:Ie,allFormatters:ce,Options:Ce,parse:Ze,gpsOnlyOptions:Ct,gps:dn,thumbnailOnlyOptions:vt,thumbnail:cn,thumbnailUrl:hn,orientationOnlyOptions:Pt,orientation:At,rotations:It,get rotateCanvas(){return ye},get rotateCss(){return we},rotation:pn});te.set("url",class extends tt{async readWhole(){this.chunked=!1;let r=await be(this.input);r instanceof ArrayBuffer?this._swapArrayBuffer(r):r instanceof Uint8Array&&this._swapBuffer(r)}async _readChunk(r,e){let t=e?r+e-1:void 0,n=this.options.httpHeaders||{};(r||t)&&(n.range=`bytes=${[r,t].join("-")}`);let i=await xt(this.input,{headers:n}),s=await i.arrayBuffer(),a=s.byteLength;if(i.status!==416)return a!==e&&(this.size=r+a),this.set(s,r,!0)}});j.prototype.getUint64=function(r){let e=this.getUint32(r),t=this.getUint32(r+4);return e<1048575?e<<32|t:typeof Ne!==void 0?(console.warn("Using BigInt because of type 64uint but JS can only handle 53b numbers."),Ne(e)<<Ne(32)|Ne(t)):void M("Trying to read 64b value but JS can only handle 53b numbers.")};class sr extends et{parseBoxes(e=0){let t=[];for(;e<this.file.byteLength-4;){let n=this.parseBoxHead(e);if(t.push(n),n.length===0)break;e+=n.length}return t}parseSubBoxes(e){e.boxes=this.parseBoxes(e.start)}findBox(e,t){return e.boxes===void 0&&this.parseSubBoxes(e),e.boxes.find(n=>n.kind===t)}parseBoxHead(e){let t=this.file.getUint32(e),n=this.file.getString(e+4,4),i=e+8;return t===1&&(t=this.file.getUint64(e+8),i+=8),{offset:e,length:t,kind:n,start:i}}parseBoxFullHead(e){if(e.version!==void 0)return;let t=this.file.getUint32(e.start);e.version=t>>24,e.start+=4}}class un extends sr{static canHandle(e,t){if(t!==0)return!1;let n=e.getUint16(2);if(n>50)return!1;let i=16,s=[];for(;i<n;)s.push(e.getString(i,4)),i+=4;return s.includes(this.type)}async parse(){let e=this.file.getUint32(0),t=this.parseBoxHead(e);for(;t.kind!=="meta";)e+=t.length,await this.file.ensureChunk(e,16),t=this.parseBoxHead(e);await this.file.ensureChunk(t.offset,t.length),this.parseBoxFullHead(t),this.parseSubBoxes(t),this.options.icc.enabled&&await this.findIcc(t),this.options.tiff.enabled&&await this.findExif(t)}async registerSegment(e,t,n){await this.file.ensureChunk(t,n);let i=this.file.subarray(t,n);this.createParser(e,i)}async findIcc(e){let t=this.findBox(e,"iprp");if(t===void 0)return;let n=this.findBox(t,"ipco");if(n===void 0)return;let i=this.findBox(n,"colr");i!==void 0&&await this.registerSegment("icc",i.offset+12,i.length)}async findExif(e){let t=this.findBox(e,"iinf");if(t===void 0)return;let n=this.findBox(e,"iloc");if(n===void 0)return;let i=this.findExifLocIdInIinf(t),s=this.findExtentInIloc(n,i);if(s===void 0)return;let[a,o]=s;await this.file.ensureChunk(a,o);let l=4+this.file.getUint32(a);a+=l,o-=l,await this.registerSegment("tiff",a,o)}findExifLocIdInIinf(e){this.parseBoxFullHead(e);let t,n,i,s,a=e.start,o=this.file.getUint16(a);for(a+=2;o--;){if(t=this.parseBoxHead(a),this.parseBoxFullHead(t),n=t.start,t.version>=2&&(i=t.version===3?4:2,s=this.file.getString(n+i+2,4),s==="Exif"))return this.file.getUintBytes(n,i);a+=t.length}}get8bits(e){let t=this.file.getUint8(e);return[t>>4,15&t]}findExtentInIloc(e,t){this.parseBoxFullHead(e);let n=e.start,[i,s]=this.get8bits(n++),[a,o]=this.get8bits(n++),l=e.version===2?4:2,h=e.version===1||e.version===2?2:0,c=o+i+s,d=e.version===2?4:2,u=this.file.getUintBytes(n,d);for(n+=d;u--;){let p=this.file.getUintBytes(n,l);n+=l+h+2+a;let _=this.file.getUint16(n);if(n+=2,p===t)return _>1&&console.warn(`ILOC box has more than one extent but we're only processing one
Please create an issue at https://github.com/MikeKovarik/exifr with this file`),[this.file.getUintBytes(n+o,i),this.file.getUintBytes(n+o+i,s)];n+=_*c}}}class fn extends un{}m(fn,"type","heic");class Gt extends un{}m(Gt,"type","avif"),ee.set("heic",fn),ee.set("avif",Gt),L(N,["ifd0","ifd1"],[[256,"ImageWidth"],[257,"ImageHeight"],[258,"BitsPerSample"],[259,"Compression"],[262,"PhotometricInterpretation"],[270,"ImageDescription"],[271,"Make"],[272,"Model"],[273,"StripOffsets"],[274,"Orientation"],[277,"SamplesPerPixel"],[278,"RowsPerStrip"],[279,"StripByteCounts"],[282,"XResolution"],[283,"YResolution"],[284,"PlanarConfiguration"],[296,"ResolutionUnit"],[301,"TransferFunction"],[305,"Software"],[306,"ModifyDate"],[315,"Artist"],[316,"HostComputer"],[317,"Predictor"],[318,"WhitePoint"],[319,"PrimaryChromaticities"],[513,"ThumbnailOffset"],[514,"ThumbnailLength"],[529,"YCbCrCoefficients"],[530,"YCbCrSubSampling"],[531,"YCbCrPositioning"],[532,"ReferenceBlackWhite"],[700,"ApplicationNotes"],[33432,"Copyright"],[33723,"IPTC"],[34665,"ExifIFD"],[34675,"ICC"],[34853,"GpsIFD"],[330,"SubIFD"],[40965,"InteropIFD"],[40091,"XPTitle"],[40092,"XPComment"],[40093,"XPAuthor"],[40094,"XPKeywords"],[40095,"XPSubject"]]),L(N,"exif",[[33434,"ExposureTime"],[33437,"FNumber"],[34850,"ExposureProgram"],[34852,"SpectralSensitivity"],[34855,"ISO"],[34858,"TimeZoneOffset"],[34859,"SelfTimerMode"],[34864,"SensitivityType"],[34865,"StandardOutputSensitivity"],[34866,"RecommendedExposureIndex"],[34867,"ISOSpeed"],[34868,"ISOSpeedLatitudeyyy"],[34869,"ISOSpeedLatitudezzz"],[36864,"ExifVersion"],[36867,"DateTimeOriginal"],[36868,"CreateDate"],[36873,"GooglePlusUploadCode"],[36880,"OffsetTime"],[36881,"OffsetTimeOriginal"],[36882,"OffsetTimeDigitized"],[37121,"ComponentsConfiguration"],[37122,"CompressedBitsPerPixel"],[37377,"ShutterSpeedValue"],[37378,"ApertureValue"],[37379,"BrightnessValue"],[37380,"ExposureCompensation"],[37381,"MaxApertureValue"],[37382,"SubjectDistance"],[37383,"MeteringMode"],[37384,"LightSource"],[37385,"Flash"],[37386,"FocalLength"],[37393,"ImageNumber"],[37394,"SecurityClassification"],[37395,"ImageHistory"],[37396,"SubjectArea"],[37500,"MakerNote"],[37510,"UserComment"],[37520,"SubSecTime"],[37521,"SubSecTimeOriginal"],[37522,"SubSecTimeDigitized"],[37888,"AmbientTemperature"],[37889,"Humidity"],[37890,"Pressure"],[37891,"WaterDepth"],[37892,"Acceleration"],[37893,"CameraElevationAngle"],[40960,"FlashpixVersion"],[40961,"ColorSpace"],[40962,"ExifImageWidth"],[40963,"ExifImageHeight"],[40964,"RelatedSoundFile"],[41483,"FlashEnergy"],[41486,"FocalPlaneXResolution"],[41487,"FocalPlaneYResolution"],[41488,"FocalPlaneResolutionUnit"],[41492,"SubjectLocation"],[41493,"ExposureIndex"],[41495,"SensingMethod"],[41728,"FileSource"],[41729,"SceneType"],[41730,"CFAPattern"],[41985,"CustomRendered"],[41986,"ExposureMode"],[41987,"WhiteBalance"],[41988,"DigitalZoomRatio"],[41989,"FocalLengthIn35mmFormat"],[41990,"SceneCaptureType"],[41991,"GainControl"],[41992,"Contrast"],[41993,"Saturation"],[41994,"Sharpness"],[41996,"SubjectDistanceRange"],[42016,"ImageUniqueID"],[42032,"OwnerName"],[42033,"SerialNumber"],[42034,"LensInfo"],[42035,"LensMake"],[42036,"LensModel"],[42037,"LensSerialNumber"],[42080,"CompositeImage"],[42081,"CompositeImageCount"],[42082,"CompositeImageExposureTimes"],[42240,"Gamma"],[59932,"Padding"],[59933,"OffsetSchema"],[65e3,"OwnerName"],[65001,"SerialNumber"],[65002,"Lens"],[65100,"RawFile"],[65101,"Converter"],[65102,"WhiteBalance"],[65105,"Exposure"],[65106,"Shadows"],[65107,"Brightness"],[65108,"Contrast"],[65109,"Saturation"],[65110,"Sharpness"],[65111,"Smoothness"],[65112,"MoireFilter"],[40965,"InteropIFD"]]),L(N,"gps",[[0,"GPSVersionID"],[1,"GPSLatitudeRef"],[2,"GPSLatitude"],[3,"GPSLongitudeRef"],[4,"GPSLongitude"],[5,"GPSAltitudeRef"],[6,"GPSAltitude"],[7,"GPSTimeStamp"],[8,"GPSSatellites"],[9,"GPSStatus"],[10,"GPSMeasureMode"],[11,"GPSDOP"],[12,"GPSSpeedRef"],[13,"GPSSpeed"],[14,"GPSTrackRef"],[15,"GPSTrack"],[16,"GPSImgDirectionRef"],[17,"GPSImgDirection"],[18,"GPSMapDatum"],[19,"GPSDestLatitudeRef"],[20,"GPSDestLatitude"],[21,"GPSDestLongitudeRef"],[22,"GPSDestLongitude"],[23,"GPSDestBearingRef"],[24,"GPSDestBearing"],[25,"GPSDestDistanceRef"],[26,"GPSDestDistance"],[27,"GPSProcessingMethod"],[28,"GPSAreaInformation"],[29,"GPSDateStamp"],[30,"GPSDifferential"],[31,"GPSHPositioningError"]]),L($,["ifd0","ifd1"],[[274,{1:"Horizontal (normal)",2:"Mirror horizontal",3:"Rotate 180",4:"Mirror vertical",5:"Mirror horizontal and rotate 270 CW",6:"Rotate 90 CW",7:"Mirror horizontal and rotate 90 CW",8:"Rotate 270 CW"}],[296,{1:"None",2:"inches",3:"cm"}]]);let De=L($,"exif",[[34850,{0:"Not defined",1:"Manual",2:"Normal program",3:"Aperture priority",4:"Shutter priority",5:"Creative program",6:"Action program",7:"Portrait mode",8:"Landscape mode"}],[37121,{0:"-",1:"Y",2:"Cb",3:"Cr",4:"R",5:"G",6:"B"}],[37383,{0:"Unknown",1:"Average",2:"CenterWeightedAverage",3:"Spot",4:"MultiSpot",5:"Pattern",6:"Partial",255:"Other"}],[37384,{0:"Unknown",1:"Daylight",2:"Fluorescent",3:"Tungsten (incandescent light)",4:"Flash",9:"Fine weather",10:"Cloudy weather",11:"Shade",12:"Daylight fluorescent (D 5700 - 7100K)",13:"Day white fluorescent (N 4600 - 5400K)",14:"Cool white fluorescent (W 3900 - 4500K)",15:"White fluorescent (WW 3200 - 3700K)",17:"Standard light A",18:"Standard light B",19:"Standard light C",20:"D55",21:"D65",22:"D75",23:"D50",24:"ISO studio tungsten",255:"Other"}],[37385,{0:"Flash did not fire",1:"Flash fired",5:"Strobe return light not detected",7:"Strobe return light detected",9:"Flash fired, compulsory flash mode",13:"Flash fired, compulsory flash mode, return light not detected",15:"Flash fired, compulsory flash mode, return light detected",16:"Flash did not fire, compulsory flash mode",24:"Flash did not fire, auto mode",25:"Flash fired, auto mode",29:"Flash fired, auto mode, return light not detected",31:"Flash fired, auto mode, return light detected",32:"No flash function",65:"Flash fired, red-eye reduction mode",69:"Flash fired, red-eye reduction mode, return light not detected",71:"Flash fired, red-eye reduction mode, return light detected",73:"Flash fired, compulsory flash mode, red-eye reduction mode",77:"Flash fired, compulsory flash mode, red-eye reduction mode, return light not detected",79:"Flash fired, compulsory flash mode, red-eye reduction mode, return light detected",89:"Flash fired, auto mode, red-eye reduction mode",93:"Flash fired, auto mode, return light not detected, red-eye reduction mode",95:"Flash fired, auto mode, return light detected, red-eye reduction mode"}],[41495,{1:"Not defined",2:"One-chip color area sensor",3:"Two-chip color area sensor",4:"Three-chip color area sensor",5:"Color sequential area sensor",7:"Trilinear sensor",8:"Color sequential linear sensor"}],[41728,{1:"Film Scanner",2:"Reflection Print Scanner",3:"Digital Camera"}],[41729,{1:"Directly photographed"}],[41985,{0:"Normal",1:"Custom",2:"HDR (no original saved)",3:"HDR (original saved)",4:"Original (for HDR)",6:"Panorama",7:"Portrait HDR",8:"Portrait"}],[41986,{0:"Auto",1:"Manual",2:"Auto bracket"}],[41987,{0:"Auto",1:"Manual"}],[41990,{0:"Standard",1:"Landscape",2:"Portrait",3:"Night",4:"Other"}],[41991,{0:"None",1:"Low gain up",2:"High gain up",3:"Low gain down",4:"High gain down"}],[41996,{0:"Unknown",1:"Macro",2:"Close",3:"Distant"}],[42080,{0:"Unknown",1:"Not a Composite Image",2:"General Composite Image",3:"Composite Image Captured While Shooting"}]]);const Ht={1:"No absolute unit of measurement",2:"Inch",3:"Centimeter"};De.set(37392,Ht),De.set(41488,Ht);const ht={0:"Normal",1:"Low",2:"High"};function jt(r){return typeof r=="object"&&r.length!==void 0?r[0]:r}function Wt(r){let e=Array.from(r).slice(1);return e[1]>15&&(e=e.map(t=>String.fromCharCode(t))),e[2]!=="0"&&e[2]!==0||e.pop(),e.join(".")}function pt(r){if(typeof r=="string"){var[e,t,n,i,s,a]=r.trim().split(/[-: ]/g).map(Number),o=new Date(e,t-1,n);return Number.isNaN(i)||Number.isNaN(s)||Number.isNaN(a)||(o.setHours(i),o.setMinutes(s),o.setSeconds(a)),Number.isNaN(+o)?r:o}}function Oe(r){if(typeof r=="string")return r;let e=[];if(r[1]===0&&r[r.length-1]===0)for(let t=0;t<r.length;t+=2)e.push(Xt(r[t+1],r[t]));else for(let t=0;t<r.length;t+=2)e.push(Xt(r[t],r[t+1]));return _e(String.fromCodePoint(...e))}function Xt(r,e){return r<<8|e}De.set(41992,ht),De.set(41993,ht),De.set(41994,ht),L(he,["ifd0","ifd1"],[[50827,function(r){return typeof r!="string"?on(r):r}],[306,pt],[40091,Oe],[40092,Oe],[40093,Oe],[40094,Oe],[40095,Oe]]),L(he,"exif",[[40960,Wt],[36864,Wt],[36867,pt],[36868,pt],[40962,jt],[40963,jt]]),L(he,"gps",[[0,r=>Array.from(r).join(".")],[7,r=>Array.from(r).join(":")]]);class ut extends Y{static canHandle(e,t){return e.getUint8(t+1)===225&&e.getUint32(t+4)===1752462448&&e.getString(t+4,20)==="http://ns.adobe.com/"}static headerLength(e,t){return e.getString(t+4,34)==="http://ns.adobe.com/xmp/extension/"?79:33}static findPosition(e,t){let n=super.findPosition(e,t);return n.multiSegment=n.extended=n.headerLength===79,n.multiSegment?(n.chunkCount=e.getUint8(t+72),n.chunkNumber=e.getUint8(t+76),e.getUint8(t+77)!==0&&n.chunkNumber++):(n.chunkCount=1/0,n.chunkNumber=-1),n}static handleMultiSegments(e){return e.map(t=>t.chunk.getString()).join("")}normalizeInput(e){return typeof e=="string"?e:j.from(e).getString()}parse(e=this.chunk){if(!this.localOptions.parse)return e;e=function(s){let a={},o={};for(let l of yn)a[l]=[],o[l]=0;return s.replace(dr,(l,h,c)=>{if(h==="<"){let d=++o[c];return a[c].push(d),`${l}#${d}`}return`${l}#${a[c].pop()}`})}(e);let t=Pe.findAll(e,"rdf","Description");t.length===0&&t.push(new Pe("rdf","Description",void 0,e));let n,i={};for(let s of t)for(let a of s.properties)n=lr(a.ns,i),gn(a,n);return function(s){let a;for(let o in s)a=s[o]=Xe(s[o]),a===void 0&&delete s[o];return Xe(s)}(i)}assignToOutput(e,t){if(this.localOptions.parse)for(let[n,i]of Object.entries(t))switch(n){case"tiff":this.assignObjectToOutput(e,"ifd0",i);break;case"exif":this.assignObjectToOutput(e,"exif",i);break;case"xmlns":break;default:this.assignObjectToOutput(e,n,i)}else e.xmp=t}}m(ut,"type","xmp"),m(ut,"multiSegment",!0),F.set("xmp",ut);class qe{static findAll(e){return mn(e,/([a-zA-Z0-9-]+):([a-zA-Z0-9-]+)=("[^"]*"|'[^']*')/gm).map(qe.unpackMatch)}static unpackMatch(e){let t=e[1],n=e[2],i=e[3].slice(1,-1);return i=_n(i),new qe(t,n,i)}constructor(e,t,n){this.ns=e,this.name=t,this.value=n}serialize(){return this.value}}class Pe{static findAll(e,t,n){if(t!==void 0||n!==void 0){t=t||"[\\w\\d-]+",n=n||"[\\w\\d-]+";var i=new RegExp(`<(${t}):(${n})(#\\d+)?((\\s+?[\\w\\d-:]+=("[^"]*"|'[^']*'))*\\s*)(\\/>|>([\\s\\S]*?)<\\/\\1:\\2\\3>)`,"gm")}else i=/<([\w\d-]+):([\w\d-]+)(#\d+)?((\s+?[\w\d-:]+=("[^"]*"|'[^']*'))*\s*)(\/>|>([\s\S]*?)<\/\1:\2\3>)/gm;return mn(e,i).map(Pe.unpackMatch)}static unpackMatch(e){let t=e[1],n=e[2],i=e[4],s=e[8];return new Pe(t,n,i,s)}constructor(e,t,n,i){this.ns=e,this.name=t,this.attrString=n,this.innerXml=i,this.attrs=qe.findAll(n),this.children=Pe.findAll(i),this.value=this.children.length===0?_n(i):void 0,this.properties=[...this.attrs,...this.children]}get isPrimitive(){return this.value!==void 0&&this.attrs.length===0&&this.children.length===0}get isListContainer(){return this.children.length===1&&this.children[0].isList}get isList(){let{ns:e,name:t}=this;return e==="rdf"&&(t==="Seq"||t==="Bag"||t==="Alt")}get isListItem(){return this.ns==="rdf"&&this.name==="li"}serialize(){if(this.properties.length===0&&this.value===void 0)return;if(this.isPrimitive)return this.value;if(this.isListContainer)return this.children[0].serialize();if(this.isList)return or(this.children.map(ar));if(this.isListItem&&this.children.length===1&&this.attrs.length===0)return this.children[0].serialize();let e={};for(let t of this.properties)gn(t,e);return this.value!==void 0&&(e.value=this.value),Xe(e)}}function gn(r,e){let t=r.serialize();t!==void 0&&(e[r.name]=t)}var ar=r=>r.serialize(),or=r=>r.length===1?r[0]:r,lr=(r,e)=>e[r]?e[r]:e[r]={};function mn(r,e){let t,n=[];if(!r)return n;for(;(t=e.exec(r))!==null;)n.push(t);return n}function _n(r){if(function(n){return n==null||n==="null"||n==="undefined"||n===""||n.trim()===""}(r))return;let e=Number(r);if(!Number.isNaN(e))return e;let t=r.toLowerCase();return t==="true"||t!=="false"&&r.trim()}const yn=["rdf:li","rdf:Seq","rdf:Bag","rdf:Alt","rdf:Description"],dr=new RegExp(`(<|\\/)(${yn.join("|")})`,"g");var wn=Object.freeze({__proto__:null,default:ir,Exifr:pe,fileParsers:ee,segmentParsers:F,fileReaders:te,tagKeys:N,tagValues:$,tagRevivers:he,createDictionary:L,extendDictionary:Se,fetchUrlAsArrayBuffer:be,readBlobAsArrayBuffer:xe,chunkedProps:oe,otherSegments:Ae,segments:ke,tiffBlocks:I,segmentsAndBlocks:le,tiffExtractables:de,inheritables:Ie,allFormatters:ce,Options:Ce,parse:Ze,gpsOnlyOptions:Ct,gps:dn,thumbnailOnlyOptions:vt,thumbnail:cn,thumbnailUrl:hn,orientationOnlyOptions:Pt,orientation:At,rotations:It,get rotateCanvas(){return ye},get rotateCss(){return we},rotation:pn});let Kt=We("fs",r=>r.promises);te.set("fs",class extends tt{async readWhole(){this.chunked=!1,this.fs=await Kt;let r=await this.fs.readFile(this.input);this._swapBuffer(r)}async readChunked(){this.chunked=!0,this.fs=await Kt,await this.open(),await this.readChunk(0,this.options.firstChunkSize)}async open(){this.fh===void 0&&(this.fh=await this.fs.open(this.input,"r"),this.size=(await this.fh.stat(this.input)).size)}async _readChunk(r,e){this.fh===void 0&&await this.open(),r+e>this.size&&(e=this.size-r);var t=this.subarray(r,e,!0);return await this.fh.read(t.dataView,0,e,r),t}async close(){if(this.fh){let r=this.fh;this.fh=void 0,await r.close()}}});te.set("base64",class extends tt{constructor(...r){super(...r),this.input=this.input.replace(/^data:([^;]+);base64,/gim,""),this.size=this.input.length/4*3,this.input.endsWith("==")?this.size-=2:this.input.endsWith("=")&&(this.size-=1)}async _readChunk(r,e){let t,n,i=this.input;r===void 0?(r=0,t=0,n=0):(t=4*Math.floor(r/3),n=r-t/4*3),e===void 0&&(e=this.size);let s=r+e,a=t+4*Math.ceil(s/3);i=i.slice(t,a);let o=Math.min(e,this.size-r);if(Qe){let l=Je.from(i,"base64").slice(n,n+o);return this.set(l,r,!0)}{let l=this.subarray(r,o,!0),h=atob(i),c=l.toUint8();for(let d=0;d<o;d++)c[d]=h.charCodeAt(n+d);return l}}});class qt extends et{static canHandle(e,t){return t===18761||t===19789}extendOptions(e){let{ifd0:t,xmp:n,iptc:i,icc:s}=e;n.enabled&&t.deps.add(700),i.enabled&&t.deps.add(33723),s.enabled&&t.deps.add(34675),t.finalizeFilters()}async parse(){let{tiff:e,xmp:t,iptc:n,icc:i}=this.options;if(e.enabled||t.enabled||n.enabled||i.enabled){let s=Math.max(gt(this.options),this.options.chunkSize);await this.file.ensureChunk(0,s),this.createParser("tiff",this.file),this.parsers.tiff.parseHeader(),await this.parsers.tiff.parseIfd0Block(),this.adaptTiffPropAsSegment("xmp"),this.adaptTiffPropAsSegment("iptc"),this.adaptTiffPropAsSegment("icc")}}adaptTiffPropAsSegment(e){if(this.parsers.tiff[e]){let t=this.parsers.tiff[e];this.injectSegment(e,t)}}}m(qt,"type","tiff"),ee.set("tiff",qt);let cr=We("zlib");const hr=["ihdr","iccp","text","itxt","exif"];class Yt extends et{constructor(...e){super(...e),m(this,"catchError",t=>this.errors.push(t)),m(this,"metaChunks",[]),m(this,"unknownChunks",[])}static canHandle(e,t){return t===35152&&e.getUint32(0)===2303741511&&e.getUint32(4)===218765834}async parse(){let{file:e}=this;await this.findPngChunksInRange(8,e.byteLength),await this.readSegments(this.metaChunks),this.findIhdr(),this.parseTextChunks(),await this.findExif().catch(this.catchError),await this.findXmp().catch(this.catchError),await this.findIcc().catch(this.catchError)}async findPngChunksInRange(e,t){let{file:n}=this;for(;e<t;){let i=n.getUint32(e),s=n.getUint32(e+4),a=n.getString(e+4,4).toLowerCase(),o=i+4+4+4,l={type:a,offset:e,length:o,start:e+4+4,size:i,marker:s};hr.includes(a)?this.metaChunks.push(l):this.unknownChunks.push(l),e+=o}}parseTextChunks(){let e=this.metaChunks.filter(t=>t.type==="text");for(let t of e){let[n,i]=this.file.getString(t.start,t.size).split("\0");this.injectKeyValToIhdr(n,i)}}injectKeyValToIhdr(e,t){let n=this.parsers.ihdr;n&&n.raw.set(e,t)}findIhdr(){let e=this.metaChunks.find(t=>t.type==="ihdr");e&&this.options.ihdr.enabled!==!1&&this.createParser("ihdr",e.chunk)}async findExif(){let e=this.metaChunks.find(t=>t.type==="exif");e&&this.injectSegment("tiff",e.chunk)}async findXmp(){let e=this.metaChunks.filter(t=>t.type==="itxt");for(let t of e)t.chunk.getString(0,17)==="XML:com.adobe.xmp"&&this.injectSegment("xmp",t.chunk)}async findIcc(){let e=this.metaChunks.find(o=>o.type==="iccp");if(!e)return;let{chunk:t}=e,n=t.getUint8Array(0,81),i=0;for(;i<80&&n[i]!==0;)i++;let s=i+2,a=t.getString(0,i);if(this.injectKeyValToIhdr("ProfileName",a),je){let o=await cr,l=t.getUint8Array(s);l=o.inflateSync(l),this.injectSegment("icc",l)}}}m(Yt,"type","png"),ee.set("png",Yt),L(N,"interop",[[1,"InteropIndex"],[2,"InteropVersion"],[4096,"RelatedImageFileFormat"],[4097,"RelatedImageWidth"],[4098,"RelatedImageHeight"]]),Se(N,"ifd0",[[11,"ProcessingSoftware"],[254,"SubfileType"],[255,"OldSubfileType"],[263,"Thresholding"],[264,"CellWidth"],[265,"CellLength"],[266,"FillOrder"],[269,"DocumentName"],[280,"MinSampleValue"],[281,"MaxSampleValue"],[285,"PageName"],[286,"XPosition"],[287,"YPosition"],[290,"GrayResponseUnit"],[297,"PageNumber"],[321,"HalftoneHints"],[322,"TileWidth"],[323,"TileLength"],[332,"InkSet"],[337,"TargetPrinter"],[18246,"Rating"],[18249,"RatingPercent"],[33550,"PixelScale"],[34264,"ModelTransform"],[34377,"PhotoshopSettings"],[50706,"DNGVersion"],[50707,"DNGBackwardVersion"],[50708,"UniqueCameraModel"],[50709,"LocalizedCameraModel"],[50736,"DNGLensInfo"],[50739,"ShadowScale"],[50740,"DNGPrivateData"],[33920,"IntergraphMatrix"],[33922,"ModelTiePoint"],[34118,"SEMInfo"],[34735,"GeoTiffDirectory"],[34736,"GeoTiffDoubleParams"],[34737,"GeoTiffAsciiParams"],[50341,"PrintIM"],[50721,"ColorMatrix1"],[50722,"ColorMatrix2"],[50723,"CameraCalibration1"],[50724,"CameraCalibration2"],[50725,"ReductionMatrix1"],[50726,"ReductionMatrix2"],[50727,"AnalogBalance"],[50728,"AsShotNeutral"],[50729,"AsShotWhiteXY"],[50730,"BaselineExposure"],[50731,"BaselineNoise"],[50732,"BaselineSharpness"],[50734,"LinearResponseLimit"],[50735,"CameraSerialNumber"],[50741,"MakerNoteSafety"],[50778,"CalibrationIlluminant1"],[50779,"CalibrationIlluminant2"],[50781,"RawDataUniqueID"],[50827,"OriginalRawFileName"],[50828,"OriginalRawFileData"],[50831,"AsShotICCProfile"],[50832,"AsShotPreProfileMatrix"],[50833,"CurrentICCProfile"],[50834,"CurrentPreProfileMatrix"],[50879,"ColorimetricReference"],[50885,"SRawType"],[50898,"PanasonicTitle"],[50899,"PanasonicTitle2"],[50931,"CameraCalibrationSig"],[50932,"ProfileCalibrationSig"],[50933,"ProfileIFD"],[50934,"AsShotProfileName"],[50936,"ProfileName"],[50937,"ProfileHueSatMapDims"],[50938,"ProfileHueSatMapData1"],[50939,"ProfileHueSatMapData2"],[50940,"ProfileToneCurve"],[50941,"ProfileEmbedPolicy"],[50942,"ProfileCopyright"],[50964,"ForwardMatrix1"],[50965,"ForwardMatrix2"],[50966,"PreviewApplicationName"],[50967,"PreviewApplicationVersion"],[50968,"PreviewSettingsName"],[50969,"PreviewSettingsDigest"],[50970,"PreviewColorSpace"],[50971,"PreviewDateTime"],[50972,"RawImageDigest"],[50973,"OriginalRawFileDigest"],[50981,"ProfileLookTableDims"],[50982,"ProfileLookTableData"],[51043,"TimeCodes"],[51044,"FrameRate"],[51058,"TStop"],[51081,"ReelName"],[51089,"OriginalDefaultFinalSize"],[51090,"OriginalBestQualitySize"],[51091,"OriginalDefaultCropSize"],[51105,"CameraLabel"],[51107,"ProfileHueSatMapEncoding"],[51108,"ProfileLookTableEncoding"],[51109,"BaselineExposureOffset"],[51110,"DefaultBlackRender"],[51111,"NewRawImageDigest"],[51112,"RawToPreviewGain"]]);let $t=[[273,"StripOffsets"],[279,"StripByteCounts"],[288,"FreeOffsets"],[289,"FreeByteCounts"],[291,"GrayResponseCurve"],[292,"T4Options"],[293,"T6Options"],[300,"ColorResponseUnit"],[320,"ColorMap"],[324,"TileOffsets"],[325,"TileByteCounts"],[326,"BadFaxLines"],[327,"CleanFaxData"],[328,"ConsecutiveBadFaxLines"],[330,"SubIFD"],[333,"InkNames"],[334,"NumberofInks"],[336,"DotRange"],[338,"ExtraSamples"],[339,"SampleFormat"],[340,"SMinSampleValue"],[341,"SMaxSampleValue"],[342,"TransferRange"],[343,"ClipPath"],[344,"XClipPathUnits"],[345,"YClipPathUnits"],[346,"Indexed"],[347,"JPEGTables"],[351,"OPIProxy"],[400,"GlobalParametersIFD"],[401,"ProfileType"],[402,"FaxProfile"],[403,"CodingMethods"],[404,"VersionYear"],[405,"ModeNumber"],[433,"Decode"],[434,"DefaultImageColor"],[435,"T82Options"],[437,"JPEGTables"],[512,"JPEGProc"],[515,"JPEGRestartInterval"],[517,"JPEGLosslessPredictors"],[518,"JPEGPointTransforms"],[519,"JPEGQTables"],[520,"JPEGDCTables"],[521,"JPEGACTables"],[559,"StripRowCounts"],[999,"USPTOMiscellaneous"],[18247,"XP_DIP_XML"],[18248,"StitchInfo"],[28672,"SonyRawFileType"],[28688,"SonyToneCurve"],[28721,"VignettingCorrection"],[28722,"VignettingCorrParams"],[28724,"ChromaticAberrationCorrection"],[28725,"ChromaticAberrationCorrParams"],[28726,"DistortionCorrection"],[28727,"DistortionCorrParams"],[29895,"SonyCropTopLeft"],[29896,"SonyCropSize"],[32781,"ImageID"],[32931,"WangTag1"],[32932,"WangAnnotation"],[32933,"WangTag3"],[32934,"WangTag4"],[32953,"ImageReferencePoints"],[32954,"RegionXformTackPoint"],[32955,"WarpQuadrilateral"],[32956,"AffineTransformMat"],[32995,"Matteing"],[32996,"DataType"],[32997,"ImageDepth"],[32998,"TileDepth"],[33300,"ImageFullWidth"],[33301,"ImageFullHeight"],[33302,"TextureFormat"],[33303,"WrapModes"],[33304,"FovCot"],[33305,"MatrixWorldToScreen"],[33306,"MatrixWorldToCamera"],[33405,"Model2"],[33421,"CFARepeatPatternDim"],[33422,"CFAPattern2"],[33423,"BatteryLevel"],[33424,"KodakIFD"],[33445,"MDFileTag"],[33446,"MDScalePixel"],[33447,"MDColorTable"],[33448,"MDLabName"],[33449,"MDSampleInfo"],[33450,"MDPrepDate"],[33451,"MDPrepTime"],[33452,"MDFileUnits"],[33589,"AdventScale"],[33590,"AdventRevision"],[33628,"UIC1Tag"],[33629,"UIC2Tag"],[33630,"UIC3Tag"],[33631,"UIC4Tag"],[33918,"IntergraphPacketData"],[33919,"IntergraphFlagRegisters"],[33921,"INGRReserved"],[34016,"Site"],[34017,"ColorSequence"],[34018,"IT8Header"],[34019,"RasterPadding"],[34020,"BitsPerRunLength"],[34021,"BitsPerExtendedRunLength"],[34022,"ColorTable"],[34023,"ImageColorIndicator"],[34024,"BackgroundColorIndicator"],[34025,"ImageColorValue"],[34026,"BackgroundColorValue"],[34027,"PixelIntensityRange"],[34028,"TransparencyIndicator"],[34029,"ColorCharacterization"],[34030,"HCUsage"],[34031,"TrapIndicator"],[34032,"CMYKEquivalent"],[34152,"AFCP_IPTC"],[34232,"PixelMagicJBIGOptions"],[34263,"JPLCartoIFD"],[34306,"WB_GRGBLevels"],[34310,"LeafData"],[34687,"TIFF_FXExtensions"],[34688,"MultiProfiles"],[34689,"SharedData"],[34690,"T88Options"],[34732,"ImageLayer"],[34750,"JBIGOptions"],[34856,"Opto-ElectricConvFactor"],[34857,"Interlace"],[34908,"FaxRecvParams"],[34909,"FaxSubAddress"],[34910,"FaxRecvTime"],[34929,"FedexEDR"],[34954,"LeafSubIFD"],[37387,"FlashEnergy"],[37388,"SpatialFrequencyResponse"],[37389,"Noise"],[37390,"FocalPlaneXResolution"],[37391,"FocalPlaneYResolution"],[37392,"FocalPlaneResolutionUnit"],[37397,"ExposureIndex"],[37398,"TIFF-EPStandardID"],[37399,"SensingMethod"],[37434,"CIP3DataFile"],[37435,"CIP3Sheet"],[37436,"CIP3Side"],[37439,"StoNits"],[37679,"MSDocumentText"],[37680,"MSPropertySetStorage"],[37681,"MSDocumentTextPosition"],[37724,"ImageSourceData"],[40965,"InteropIFD"],[40976,"SamsungRawPointersOffset"],[40977,"SamsungRawPointersLength"],[41217,"SamsungRawByteOrder"],[41218,"SamsungRawUnknown"],[41484,"SpatialFrequencyResponse"],[41485,"Noise"],[41489,"ImageNumber"],[41490,"SecurityClassification"],[41491,"ImageHistory"],[41494,"TIFF-EPStandardID"],[41995,"DeviceSettingDescription"],[42112,"GDALMetadata"],[42113,"GDALNoData"],[44992,"ExpandSoftware"],[44993,"ExpandLens"],[44994,"ExpandFilm"],[44995,"ExpandFilterLens"],[44996,"ExpandScanner"],[44997,"ExpandFlashLamp"],[46275,"HasselbladRawImage"],[48129,"PixelFormat"],[48130,"Transformation"],[48131,"Uncompressed"],[48132,"ImageType"],[48256,"ImageWidth"],[48257,"ImageHeight"],[48258,"WidthResolution"],[48259,"HeightResolution"],[48320,"ImageOffset"],[48321,"ImageByteCount"],[48322,"AlphaOffset"],[48323,"AlphaByteCount"],[48324,"ImageDataDiscard"],[48325,"AlphaDataDiscard"],[50215,"OceScanjobDesc"],[50216,"OceApplicationSelector"],[50217,"OceIDNumber"],[50218,"OceImageLogic"],[50255,"Annotations"],[50459,"HasselbladExif"],[50547,"OriginalFileName"],[50560,"USPTOOriginalContentType"],[50656,"CR2CFAPattern"],[50710,"CFAPlaneColor"],[50711,"CFALayout"],[50712,"LinearizationTable"],[50713,"BlackLevelRepeatDim"],[50714,"BlackLevel"],[50715,"BlackLevelDeltaH"],[50716,"BlackLevelDeltaV"],[50717,"WhiteLevel"],[50718,"DefaultScale"],[50719,"DefaultCropOrigin"],[50720,"DefaultCropSize"],[50733,"BayerGreenSplit"],[50737,"ChromaBlurRadius"],[50738,"AntiAliasStrength"],[50752,"RawImageSegmentation"],[50780,"BestQualityScale"],[50784,"AliasLayerMetadata"],[50829,"ActiveArea"],[50830,"MaskedAreas"],[50935,"NoiseReductionApplied"],[50974,"SubTileBlockSize"],[50975,"RowInterleaveFactor"],[51008,"OpcodeList1"],[51009,"OpcodeList2"],[51022,"OpcodeList3"],[51041,"NoiseProfile"],[51114,"CacheVersion"],[51125,"DefaultUserCrop"],[51157,"NikonNEFInfo"],[65024,"KdcIFD"]];Se(N,"ifd0",$t),Se(N,"exif",$t),L($,"gps",[[23,{M:"Magnetic North",T:"True North"}],[25,{K:"Kilometers",M:"Miles",N:"Nautical Miles"}]]);class ft extends Y{static canHandle(e,t){return e.getUint8(t+1)===224&&e.getUint32(t+4)===1246120262&&e.getUint8(t+8)===0}parse(){return this.parseTags(),this.translate(),this.output}parseTags(){this.raw=new Map([[0,this.chunk.getUint16(0)],[2,this.chunk.getUint8(2)],[3,this.chunk.getUint16(3)],[5,this.chunk.getUint16(5)],[7,this.chunk.getUint8(7)],[8,this.chunk.getUint8(8)]])}}m(ft,"type","jfif"),m(ft,"headerLength",9),F.set("jfif",ft),L(N,"jfif",[[0,"JFIFVersion"],[2,"ResolutionUnit"],[3,"XResolution"],[5,"YResolution"],[7,"ThumbnailWidth"],[8,"ThumbnailHeight"]]);class Jt extends Y{parse(){return this.parseTags(),this.translate(),this.output}parseTags(){this.raw=new Map([[0,this.chunk.getUint32(0)],[4,this.chunk.getUint32(4)],[8,this.chunk.getUint8(8)],[9,this.chunk.getUint8(9)],[10,this.chunk.getUint8(10)],[11,this.chunk.getUint8(11)],[12,this.chunk.getUint8(12)],...Array.from(this.raw)])}}m(Jt,"type","ihdr"),F.set("ihdr",Jt),L(N,"ihdr",[[0,"ImageWidth"],[4,"ImageHeight"],[8,"BitDepth"],[9,"ColorType"],[10,"Compression"],[11,"Filter"],[12,"Interlace"]]),L($,"ihdr",[[9,{0:"Grayscale",2:"RGB",3:"Palette",4:"Grayscale with Alpha",6:"RGB with Alpha",DEFAULT:"Unknown"}],[10,{0:"Deflate/Inflate",DEFAULT:"Unknown"}],[11,{0:"Adaptive",DEFAULT:"Unknown"}],[12,{0:"Noninterlaced",1:"Adam7 Interlace",DEFAULT:"Unknown"}]]);class ze extends Y{static canHandle(e,t){return e.getUint8(t+1)===226&&e.getUint32(t+4)===1229144927}static findPosition(e,t){let n=super.findPosition(e,t);return n.chunkNumber=e.getUint8(t+16),n.chunkCount=e.getUint8(t+17),n.multiSegment=n.chunkCount>1,n}static handleMultiSegments(e){return function(t){let n=function(i){let s=i[0].constructor,a=0;for(let h of i)a+=h.length;let o=new s(a),l=0;for(let h of i)o.set(h,l),l+=h.length;return o}(t.map(i=>i.chunk.toUint8()));return new j(n)}(e)}parse(){return this.raw=new Map,this.parseHeader(),this.parseTags(),this.translate(),this.output}parseHeader(){let{raw:e}=this;this.chunk.byteLength<84&&M("ICC header is too short");for(let[t,n]of Object.entries(pr)){t=parseInt(t,10);let i=n(this.chunk,t);i!=="\0\0\0\0"&&e.set(t,i)}}parseTags(){let e,t,n,i,s,{raw:a}=this,o=this.chunk.getUint32(128),l=132,h=this.chunk.byteLength;for(;o--;){if(e=this.chunk.getString(l,4),t=this.chunk.getUint32(l+4),n=this.chunk.getUint32(l+8),i=this.chunk.getString(t,4),t+n>h)return void console.warn("reached the end of the first ICC chunk. Enable options.tiff.multiSegment to read all ICC segments.");s=this.parseTag(i,t,n),s!==void 0&&s!=="\0\0\0\0"&&a.set(e,s),l+=12}}parseTag(e,t,n){switch(e){case"desc":return this.parseDesc(t);case"mluc":return this.parseMluc(t);case"text":return this.parseText(t,n);case"sig ":return this.parseSig(t)}if(!(t+n>this.chunk.byteLength))return this.chunk.getUint8Array(t,n)}parseDesc(e){let t=this.chunk.getUint32(e+8)-1;return _e(this.chunk.getString(e+12,t))}parseText(e,t){return _e(this.chunk.getString(e+8,t-8))}parseSig(e){return _e(this.chunk.getString(e+8,4))}parseMluc(e){let{chunk:t}=this,n=t.getUint32(e+8),i=t.getUint32(e+12),s=e+16,a=[];for(let o=0;o<n;o++){let l=t.getString(s+0,2),h=t.getString(s+2,2),c=t.getUint32(s+4),d=t.getUint32(s+8)+e,u=_e(t.getUnicodeString(d,c));a.push({lang:l,country:h,text:u}),s+=i}return n===1?a[0].text:a}translateValue(e,t){return typeof e=="string"?t[e]||t[e.toLowerCase()]||e:t[e]||e}}m(ze,"type","icc"),m(ze,"multiSegment",!0),m(ze,"headerLength",18);const pr={4:ie,8:function(r,e){return[r.getUint8(e),r.getUint8(e+1)>>4,r.getUint8(e+1)%16].map(t=>t.toString(10)).join(".")},12:ie,16:ie,20:ie,24:function(r,e){const t=r.getUint16(e),n=r.getUint16(e+2)-1,i=r.getUint16(e+4),s=r.getUint16(e+6),a=r.getUint16(e+8),o=r.getUint16(e+10);return new Date(Date.UTC(t,n,i,s,a,o))},36:ie,40:ie,48:ie,52:ie,64:(r,e)=>r.getUint32(e),80:ie};function ie(r,e){return _e(r.getString(e,4))}F.set("icc",ze),L(N,"icc",[[4,"ProfileCMMType"],[8,"ProfileVersion"],[12,"ProfileClass"],[16,"ColorSpaceData"],[20,"ProfileConnectionSpace"],[24,"ProfileDateTime"],[36,"ProfileFileSignature"],[40,"PrimaryPlatform"],[44,"CMMFlags"],[48,"DeviceManufacturer"],[52,"DeviceModel"],[56,"DeviceAttributes"],[64,"RenderingIntent"],[68,"ConnectionSpaceIlluminant"],[80,"ProfileCreator"],[84,"ProfileID"],["Header","ProfileHeader"],["MS00","WCSProfiles"],["bTRC","BlueTRC"],["bXYZ","BlueMatrixColumn"],["bfd","UCRBG"],["bkpt","MediaBlackPoint"],["calt","CalibrationDateTime"],["chad","ChromaticAdaptation"],["chrm","Chromaticity"],["ciis","ColorimetricIntentImageState"],["clot","ColorantTableOut"],["clro","ColorantOrder"],["clrt","ColorantTable"],["cprt","ProfileCopyright"],["crdi","CRDInfo"],["desc","ProfileDescription"],["devs","DeviceSettings"],["dmdd","DeviceModelDesc"],["dmnd","DeviceMfgDesc"],["dscm","ProfileDescriptionML"],["fpce","FocalPlaneColorimetryEstimates"],["gTRC","GreenTRC"],["gXYZ","GreenMatrixColumn"],["gamt","Gamut"],["kTRC","GrayTRC"],["lumi","Luminance"],["meas","Measurement"],["meta","Metadata"],["mmod","MakeAndModel"],["ncl2","NamedColor2"],["ncol","NamedColor"],["ndin","NativeDisplayInfo"],["pre0","Preview0"],["pre1","Preview1"],["pre2","Preview2"],["ps2i","PS2RenderingIntent"],["ps2s","PostScript2CSA"],["psd0","PostScript2CRD0"],["psd1","PostScript2CRD1"],["psd2","PostScript2CRD2"],["psd3","PostScript2CRD3"],["pseq","ProfileSequenceDesc"],["psid","ProfileSequenceIdentifier"],["psvm","PS2CRDVMSize"],["rTRC","RedTRC"],["rXYZ","RedMatrixColumn"],["resp","OutputResponse"],["rhoc","ReflectionHardcopyOrigColorimetry"],["rig0","PerceptualRenderingIntentGamut"],["rig2","SaturationRenderingIntentGamut"],["rpoc","ReflectionPrintOutputColorimetry"],["sape","SceneAppearanceEstimates"],["scoe","SceneColorimetryEstimates"],["scrd","ScreeningDesc"],["scrn","Screening"],["targ","CharTarget"],["tech","Technology"],["vcgt","VideoCardGamma"],["view","ViewingConditions"],["vued","ViewingCondDesc"],["wtpt","MediaWhitePoint"]]);const Be={"4d2p":"Erdt Systems",AAMA:"Aamazing Technologies",ACER:"Acer",ACLT:"Acolyte Color Research",ACTI:"Actix Sytems",ADAR:"Adara Technology",ADBE:"Adobe",ADI:"ADI Systems",AGFA:"Agfa Graphics",ALMD:"Alps Electric",ALPS:"Alps Electric",ALWN:"Alwan Color Expertise",AMTI:"Amiable Technologies",AOC:"AOC International",APAG:"Apago",APPL:"Apple Computer",AST:"AST","AT&T":"AT&T",BAEL:"BARBIERI electronic",BRCO:"Barco NV",BRKP:"Breakpoint",BROT:"Brother",BULL:"Bull",BUS:"Bus Computer Systems","C-IT":"C-Itoh",CAMR:"Intel",CANO:"Canon",CARR:"Carroll Touch",CASI:"Casio",CBUS:"Colorbus PL",CEL:"Crossfield",CELx:"Crossfield",CGS:"CGS Publishing Technologies International",CHM:"Rochester Robotics",CIGL:"Colour Imaging Group, London",CITI:"Citizen",CL00:"Candela",CLIQ:"Color IQ",CMCO:"Chromaco",CMiX:"CHROMiX",COLO:"Colorgraphic Communications",COMP:"Compaq",COMp:"Compeq/Focus Technology",CONR:"Conrac Display Products",CORD:"Cordata Technologies",CPQ:"Compaq",CPRO:"ColorPro",CRN:"Cornerstone",CTX:"CTX International",CVIS:"ColorVision",CWC:"Fujitsu Laboratories",DARI:"Darius Technology",DATA:"Dataproducts",DCP:"Dry Creek Photo",DCRC:"Digital Contents Resource Center, Chung-Ang University",DELL:"Dell Computer",DIC:"Dainippon Ink and Chemicals",DICO:"Diconix",DIGI:"Digital","DL&C":"Digital Light & Color",DPLG:"Doppelganger",DS:"Dainippon Screen",DSOL:"DOOSOL",DUPN:"DuPont",EPSO:"Epson",ESKO:"Esko-Graphics",ETRI:"Electronics and Telecommunications Research Institute",EVER:"Everex Systems",EXAC:"ExactCODE",Eizo:"Eizo",FALC:"Falco Data Products",FF:"Fuji Photo Film",FFEI:"FujiFilm Electronic Imaging",FNRD:"Fnord Software",FORA:"Fora",FORE:"Forefront Technology",FP:"Fujitsu",FPA:"WayTech Development",FUJI:"Fujitsu",FX:"Fuji Xerox",GCC:"GCC Technologies",GGSL:"Global Graphics Software",GMB:"Gretagmacbeth",GMG:"GMG",GOLD:"GoldStar Technology",GOOG:"Google",GPRT:"Giantprint",GTMB:"Gretagmacbeth",GVC:"WayTech Development",GW2K:"Sony",HCI:"HCI",HDM:"Heidelberger Druckmaschinen",HERM:"Hermes",HITA:"Hitachi America",HP:"Hewlett-Packard",HTC:"Hitachi",HiTi:"HiTi Digital",IBM:"IBM",IDNT:"Scitex",IEC:"Hewlett-Packard",IIYA:"Iiyama North America",IKEG:"Ikegami Electronics",IMAG:"Image Systems",IMI:"Ingram Micro",INTC:"Intel",INTL:"N/A (INTL)",INTR:"Intra Electronics",IOCO:"Iocomm International Technology",IPS:"InfoPrint Solutions Company",IRIS:"Scitex",ISL:"Ichikawa Soft Laboratory",ITNL:"N/A (ITNL)",IVM:"IVM",IWAT:"Iwatsu Electric",Idnt:"Scitex",Inca:"Inca Digital Printers",Iris:"Scitex",JPEG:"Joint Photographic Experts Group",JSFT:"Jetsoft Development",JVC:"JVC Information Products",KART:"Scitex",KFC:"KFC Computek Components",KLH:"KLH Computers",KMHD:"Konica Minolta",KNCA:"Konica",KODA:"Kodak",KYOC:"Kyocera",Kart:"Scitex",LCAG:"Leica",LCCD:"Leeds Colour",LDAK:"Left Dakota",LEAD:"Leading Technology",LEXM:"Lexmark International",LINK:"Link Computer",LINO:"Linotronic",LITE:"Lite-On",Leaf:"Leaf",Lino:"Linotronic",MAGC:"Mag Computronic",MAGI:"MAG Innovision",MANN:"Mannesmann",MICN:"Micron Technology",MICR:"Microtek",MICV:"Microvitec",MINO:"Minolta",MITS:"Mitsubishi Electronics America",MITs:"Mitsuba",MNLT:"Minolta",MODG:"Modgraph",MONI:"Monitronix",MONS:"Monaco Systems",MORS:"Morse Technology",MOTI:"Motive Systems",MSFT:"Microsoft",MUTO:"MUTOH INDUSTRIES",Mits:"Mitsubishi Electric",NANA:"NANAO",NEC:"NEC",NEXP:"NexPress Solutions",NISS:"Nissei Sangyo America",NKON:"Nikon",NONE:"none",OCE:"Oce Technologies",OCEC:"OceColor",OKI:"Oki",OKID:"Okidata",OKIP:"Okidata",OLIV:"Olivetti",OLYM:"Olympus",ONYX:"Onyx Graphics",OPTI:"Optiquest",PACK:"Packard Bell",PANA:"Matsushita Electric Industrial",PANT:"Pantone",PBN:"Packard Bell",PFU:"PFU",PHIL:"Philips Consumer Electronics",PNTX:"HOYA",POne:"Phase One A/S",PREM:"Premier Computer Innovations",PRIN:"Princeton Graphic Systems",PRIP:"Princeton Publishing Labs",QLUX:"Hong Kong",QMS:"QMS",QPCD:"QPcard AB",QUAD:"QuadLaser",QUME:"Qume",RADI:"Radius",RDDx:"Integrated Color Solutions",RDG:"Roland DG",REDM:"REDMS Group",RELI:"Relisys",RGMS:"Rolf Gierling Multitools",RICO:"Ricoh",RNLD:"Edmund Ronald",ROYA:"Royal",RPC:"Ricoh Printing Systems",RTL:"Royal Information Electronics",SAMP:"Sampo",SAMS:"Samsung",SANT:"Jaime Santana Pomares",SCIT:"Scitex",SCRN:"Dainippon Screen",SDP:"Scitex",SEC:"Samsung",SEIK:"Seiko Instruments",SEIk:"Seikosha",SGUY:"ScanGuy.com",SHAR:"Sharp Laboratories",SICC:"International Color Consortium",SONY:"Sony",SPCL:"SpectraCal",STAR:"Star",STC:"Sampo Technology",Scit:"Scitex",Sdp:"Scitex",Sony:"Sony",TALO:"Talon Technology",TAND:"Tandy",TATU:"Tatung",TAXA:"TAXAN America",TDS:"Tokyo Denshi Sekei",TECO:"TECO Information Systems",TEGR:"Tegra",TEKT:"Tektronix",TI:"Texas Instruments",TMKR:"TypeMaker",TOSB:"Toshiba",TOSH:"Toshiba",TOTK:"TOTOKU ELECTRIC",TRIU:"Triumph",TSBT:"Toshiba",TTX:"TTX Computer Products",TVM:"TVM Professional Monitor",TW:"TW Casper",ULSX:"Ulead Systems",UNIS:"Unisys",UTZF:"Utz Fehlau & Sohn",VARI:"Varityper",VIEW:"Viewsonic",VISL:"Visual communication",VIVO:"Vivo Mobile Communication",WANG:"Wang",WLBR:"Wilbur Imaging",WTG2:"Ware To Go",WYSE:"WYSE Technology",XERX:"Xerox",XRIT:"X-Rite",ZRAN:"Zoran",Zebr:"Zebra Technologies",appl:"Apple Computer",bICC:"basICColor",berg:"bergdesign",ceyd:"Integrated Color Solutions",clsp:"MacDermid ColorSpan",ds:"Dainippon Screen",dupn:"DuPont",ffei:"FujiFilm Electronic Imaging",flux:"FluxData",iris:"Scitex",kart:"Scitex",lcms:"Little CMS",lino:"Linotronic",none:"none",ob4d:"Erdt Systems",obic:"Medigraph",quby:"Qubyx Sarl",scit:"Scitex",scrn:"Dainippon Screen",sdp:"Scitex",siwi:"SIWI GRAFIKA",yxym:"YxyMaster"},Qt={scnr:"Scanner",mntr:"Monitor",prtr:"Printer",link:"Device Link",abst:"Abstract",spac:"Color Space Conversion Profile",nmcl:"Named Color",cenc:"ColorEncodingSpace profile",mid:"MultiplexIdentification profile",mlnk:"MultiplexLink profile",mvis:"MultiplexVisualization profile",nkpf:"Nikon Input Device Profile (NON-STANDARD!)"};L($,"icc",[[4,Be],[12,Qt],[40,Object.assign({},Be,Qt)],[48,Be],[80,Be],[64,{0:"Perceptual",1:"Relative Colorimetric",2:"Saturation",3:"Absolute Colorimetric"}],["tech",{amd:"Active Matrix Display",crt:"Cathode Ray Tube Display",kpcd:"Photo CD",pmd:"Passive Matrix Display",dcam:"Digital Camera",dcpj:"Digital Cinema Projector",dmpc:"Digital Motion Picture Camera",dsub:"Dye Sublimation Printer",epho:"Electrophotographic Printer",esta:"Electrostatic Printer",flex:"Flexography",fprn:"Film Writer",fscn:"Film Scanner",grav:"Gravure",ijet:"Ink Jet Printer",imgs:"Photo Image Setter",mpfr:"Motion Picture Film Recorder",mpfs:"Motion Picture Film Scanner",offs:"Offset Lithography",pjtv:"Projection Television",rpho:"Photographic Paper Printer",rscn:"Reflective Scanner",silk:"Silkscreen",twax:"Thermal Wax Printer",vidc:"Video Camera",vidm:"Video Monitor"}]]);class Ve extends Y{static canHandle(e,t,n){return e.getUint8(t+1)===237&&e.getString(t+4,9)==="Photoshop"&&this.containsIptc8bim(e,t,n)!==void 0}static headerLength(e,t,n){let i,s=this.containsIptc8bim(e,t,n);if(s!==void 0)return i=e.getUint8(t+s+7),i%2!=0&&(i+=1),i===0&&(i=4),s+8+i}static containsIptc8bim(e,t,n){for(let i=0;i<n;i++)if(this.isIptcSegmentHead(e,t+i))return i}static isIptcSegmentHead(e,t){return e.getUint8(t)===56&&e.getUint32(t)===943868237&&e.getUint16(t+4)===1028}parse(){let{raw:e}=this,t=this.chunk.byteLength-1,n=!1;for(let i=0;i<t;i++)if(this.chunk.getUint8(i)===28&&this.chunk.getUint8(i+1)===2){n=!0;let s=this.chunk.getUint16(i+3),a=this.chunk.getUint8(i+2),o=this.chunk.getLatin1String(i+5,s);e.set(a,this.pluralizeValue(e.get(a),o)),i+=4+s}else if(n)break;return this.translate(),this.output}pluralizeValue(e,t){return e!==void 0?e instanceof Array?(e.push(t),e):[e,t]:t}}m(Ve,"type","iptc"),m(Ve,"translateValues",!1),m(Ve,"reviveValues",!1),F.set("iptc",Ve),L(N,"iptc",[[0,"ApplicationRecordVersion"],[3,"ObjectTypeReference"],[4,"ObjectAttributeReference"],[5,"ObjectName"],[7,"EditStatus"],[8,"EditorialUpdate"],[10,"Urgency"],[12,"SubjectReference"],[15,"Category"],[20,"SupplementalCategories"],[22,"FixtureIdentifier"],[25,"Keywords"],[26,"ContentLocationCode"],[27,"ContentLocationName"],[30,"ReleaseDate"],[35,"ReleaseTime"],[37,"ExpirationDate"],[38,"ExpirationTime"],[40,"SpecialInstructions"],[42,"ActionAdvised"],[45,"ReferenceService"],[47,"ReferenceDate"],[50,"ReferenceNumber"],[55,"DateCreated"],[60,"TimeCreated"],[62,"DigitalCreationDate"],[63,"DigitalCreationTime"],[65,"OriginatingProgram"],[70,"ProgramVersion"],[75,"ObjectCycle"],[80,"Byline"],[85,"BylineTitle"],[90,"City"],[92,"Sublocation"],[95,"State"],[100,"CountryCode"],[101,"Country"],[103,"OriginalTransmissionReference"],[105,"Headline"],[110,"Credit"],[115,"Source"],[116,"CopyrightNotice"],[118,"Contact"],[120,"Caption"],[121,"LocalCaption"],[122,"Writer"],[125,"RasterizedCaption"],[130,"ImageType"],[131,"ImageOrientation"],[135,"LanguageIdentifier"],[150,"AudioType"],[151,"AudioSamplingRate"],[152,"AudioSamplingResolution"],[153,"AudioDuration"],[154,"AudioOutcue"],[184,"JobID"],[185,"MasterDocumentID"],[186,"ShortDocumentID"],[187,"UniqueDocumentID"],[188,"OwnerID"],[200,"ObjectPreviewFileFormat"],[201,"ObjectPreviewFileVersion"],[202,"ObjectPreviewData"],[221,"Prefs"],[225,"ClassifyState"],[228,"SimilarityIndex"],[230,"DocumentNotes"],[231,"DocumentHistory"],[232,"ExifCameraInfo"],[255,"CatalogSets"]]),L($,"iptc",[[10,{0:"0 (reserved)",1:"1 (most urgent)",2:"2",3:"3",4:"4",5:"5 (normal urgency)",6:"6",7:"7",8:"8 (least urgent)",9:"9 (user-defined priority)"}],[75,{a:"Morning",b:"Both Morning and Evening",p:"Evening"}],[131,{L:"Landscape",P:"Portrait",S:"Square"}]]);function ur(r,e,t){const n=t.type===1||t.type===2||t.type===7?1:t.type===3?2:t.type===4||t.type===9||t.type===11||t.type===13?4:t.type===5||t.type===10||t.type===12?8:0;if(n===0)throw new Error(`Unsupported TIFF field type: ${t.type}`);const s=t.count*n<=4?t.valueFieldOffset:r.getUint32(t.valueFieldOffset,e),a=[];for(let o=0;o<t.count;o++){const l=s+o*n;if(l<0||l+n>r.byteLength)throw new Error("Invalid TIFF field offset");if(t.type===1)a.push(r.getUint8(l));else if(t.type===2||t.type===7)a.push(r.getUint8(l));else if(t.type===3)a.push(r.getUint16(l,e));else if(t.type===5){const h=r.getUint32(l+4,e);a.push(h?r.getUint32(l,e)/h:0)}else if(t.type===13)a.push(r.getUint32(l,e));else if(t.type===9)a.push(r.getInt32(l,e));else if(t.type===10){const h=r.getInt32(l+4,e);a.push(h?r.getInt32(l,e)/h:0)}else t.type===11?a.push(r.getFloat32(l,e)):t.type===12?a.push(r.getFloat64(l,e)):a.push(r.getUint32(l,e))}return a}function Zt(r,e,t){if(!Number.isFinite(t)||t<=0||t+2>r.byteLength)throw new Error("Invalid TIFF IFD offset");const n=r.getUint16(t,e);if(t+2+n*12+4>r.byteLength)throw new Error("Corrupt TIFF IFD");const s=new Map;for(let a=0;a<n;a++){const o=t+2+a*12,l=r.getUint16(o,e),h=r.getUint16(o+2,e),c=r.getUint32(o+4,e);s.set(l,{type:h,count:c,valueFieldOffset:o+8})}return s}function H(r,e,t,n){const i=t.get(n);return i?ur(r,e,i):[]}const bn=(...r)=>{for(const e of r){if(Array.isArray(e)||ArrayBuffer.isView(e)){const n=Array.from(e).map(Number).filter(i=>Number.isFinite(i)&&i>0);if(n.length)return Math.max(...n);continue}const t=Number(e);if(Number.isFinite(t)&&t>0)return t}return null},fr=r=>{const t=(Array.isArray(r)||ArrayBuffer.isView(r)?Array.from(r).map(Number):[Number(r)]).map(n=>Number.isFinite(n)?n:0);return t.length>=3?[t[0]||0,t[1]||0,t[1]||0,t[2]||0]:t.length===1?[t[0]||0,t[0]||0,t[0]||0,t[0]||0]:[0,0,0,0]};function gr(r,e={}){var k,W,q,X,x,B;if(r.byteLength<8)return null;const t=new DataView(r.buffer,r.byteOffset,r.byteLength),n=t.getUint16(0,!1),i=n===18761;if(!i&&n!==19789||t.getUint16(2,i)!==42)return null;const s=t.getUint32(4,i);let a;try{a=Zt(t,i,s)}catch{return null}const o=new Set([s]);for(const f of H(t,i,a,330))f>0&&o.add(f);let l=null;for(const f of o)try{const C=f===s?a:Zt(t,i,f),D=H(t,i,C,256)[0]||0,S=H(t,i,C,257)[0]||0,w=H(t,i,C,258),R=H(t,i,C,277)[0]||w.length||1,J=w.length===1?new Array(R).fill(w[0]):w,ue=H(t,i,C,259)[0]||1,ne=H(t,i,C,262)[0]||0,Q=H(t,i,C,284)[0]||1,V=H(t,i,C,339),nt=V.length===0?new Array(R).fill(1):V.length===1?new Array(R).fill(V[0]):V,ve=H(t,i,C,273),Me=H(t,i,C,279),rt=H(t,i,C,278)[0]||S,it=nt.slice(0,R).every(se=>se===0||se===1);if(!(ne===34892&&ue===1&&Q===1&&D>0&&S>0&&R>=3&&J.length>=R&&J.slice(0,R).every(se=>se===16)&&it&&ve.length>0&&ve.length===Me.length))continue;(!l||D*S>l.width*l.height)&&(l={entries:C,width:D,height:S,samplesPerPixel:R,rowsPerStrip:rt,stripOffsets:ve,stripByteCounts:Me,bitsPerSample:J,photometric:ne})}catch{continue}if(!l)return null;const{entries:h,width:c,height:d,samplesPerPixel:u,rowsPerStrip:p,stripOffsets:_,stripByteCounts:g,bitsPerSample:y,photometric:b}=l,v=c*d;if(!Number.isSafeInteger(v)||v<=0)return null;const A=new Uint16Array(v*3);for(let f=0;f<_.length;f++){const C=f*p;if(C>=d)break;const D=Math.min(p,d-C),S=D*c*u*2,w=_[f],R=g[f];if(w<0||R<S||w+S>r.byteLength)throw new Error("Invalid LinearRaw DNG strip bounds");const J=C*c*3;if(u===3&&i&&!(r.byteOffset+w&1)){A.set(new Uint16Array(r.buffer,r.byteOffset+w,D*c*3),J);continue}const ue=new DataView(r.buffer,r.byteOffset+w,S);let ne=0,Q=J;for(let V=0;V<D*c;V++)A[Q++]=ue.getUint16(ne,i),A[Q++]=ue.getUint16(ne+2,i),A[Q++]=ue.getUint16(ne+4,i),ne+=u*2}const T=H(t,i,h,50714),P=((W=(k=e==null?void 0:e.color_data)==null?void 0:k.dng_levels)==null?void 0:W.dng_cblack)||((q=e==null?void 0:e.color_data)==null?void 0:q.cblack_rawpy_style)||(e==null?void 0:e.black_level_per_channel)||(e==null?void 0:e.cblack),O=fr(T.length?T:P),E=H(t,i,h,50717),U=bn(E,(x=(X=e==null?void 0:e.color_data)==null?void 0:X.dng_levels)==null?void 0:x.dng_whitelevel,(B=e==null?void 0:e.color_data)==null?void 0:B.maximum,e==null?void 0:e.white_level)||65535;return{data:A,width:c,height:d,bayerPattern:"",blackLevels:O,whiteLevel:U,metadata:{...e,format:"DNG_LINEAR_RAW_RGB",description:"Uncompressed DNG LinearRaw RGB",linearRawDngDecoder:!0,bitsPerSample:y.slice(0,u),samplesPerPixel:u,photometric:b},isThreePlane:!0,threePlaneTransfer:"linear"}}async function mr(r){const e=gr(r);if(!e)return null;let t={};const n=await bt(),i=new n;try{await i.open(r,{}),t=await i.metadata(!0)}catch(s){console.warn("LinearRaw DNG metadata enrichment failed",s)}finally{try{i.delete?i.delete():i.close()}catch{}}try{const s=await wn.parse(r.buffer);s&&(t={...t,...s})}catch(s){console.warn("exifr parsing failed for LinearRaw DNG",s)}return e.metadata={...t,...e.metadata},e}const _r=async r=>{var i,s,a,o,l,h,c,d,u,p,_;const e=await mr(r);if(e)return e;const t=await bt(),n=new t;try{if(await n.open(r,{}),typeof n.getRawImage!="function")throw new Error("WASM mismatch");let g={};try{g=await n.metadata(!0)}catch(x){console.warn("Metadata error before raw extraction",x)}const y=n.getRawImage(),b=y.data instanceof Uint16Array?y.data:new Uint16Array(y.data);let v={...g};try{const x=await wn.parse(r.buffer);x&&(v={...v,...x})}catch(x){console.warn("exifr parsing failed for RAW buffer",x)}const A=g.filters||((i=g.idata)==null?void 0:i.filters)||0,T=g.colors||((s=g.idata)==null?void 0:s.colors)||0,P=A===0&&T===3,O=A===9;let E=[0,0,0,0],U=null,k;if(n.getBlackLevels)try{const x=n.getBlackLevels();U=x,k=Ft(x)||Ft((a=g.color_data)==null?void 0:a.black_level_model)||void 0;const B=me(k==null?void 0:k.siteBaseLevels);B&&(E=B)}catch(x){console.warn("getBlackLevels binding failed",x)}if(!k){const x=Number((U==null?void 0:U.black)??((o=g.color_data)==null?void 0:o.black)??0)||0,B=me((U==null?void 0:U.cblack)||((l=g.color_data)==null?void 0:l.cblack_rawpy_style)||g.black_level_per_channel||g.cblack||((h=g.color)==null?void 0:h.cblack));if(B){const f=String(y.bayerPattern||g.cfa_pattern||"RGGB").toUpperCase();let C=0;E=[0,1,2,3].map(D=>{const S=f[D],w=S==="R"?0:S==="B"?2:++C===1?1:3;return Math.max(0,x+B[w])})}else E=[x,x,x,x]}const W=Number.isFinite(Number(y.bits))&&Number(y.bits)>0?Math.pow(2,Number(y.bits))-1:null,q=bn(g.white_level,(d=(c=g.color_data)==null?void 0:c.dng_levels)==null?void 0:d.dng_whitelevel,(u=g.color_data)==null?void 0:u.maximum,U==null?void 0:U.maximum,(p=g.color_data)==null?void 0:p.fmaximum,(_=g.color_data)==null?void 0:_.data_maximum,W)||16383,X={data:b,width:y.width,height:y.height,bayerPattern:y.bayerPattern||"",blackLevels:E,blackLevelModel:k,whiteLevel:q,metadata:v,isThreePlane:P,threePlaneTransfer:P?"linear":void 0,isXTrans:O};return Gn(X,zn(X,v)),X}finally{n.delete?n.delete():n.close()}};async function yr(r){if(On(r)){const n=await Un(r);if(!n)throw new Error("Sony cRAW HQ decoder did not return image data.");return n.rawImageData}const t=new Uint8Array(r);return _r(t)}function wr(r,e,t,n,i,s=!1,a=!1,o={rows:!1,rowGroup:0,cols:!1,colGroup:0},l=!0,h=1,c=!1){const d=Math.max(0,Math.floor(n.x)),u=Math.max(0,Math.floor(n.y)),p=Math.min(e,Math.ceil(n.x+n.w)),_=Math.min(t,Math.ceil(n.y+n.h)),g=p-d,y=_-u;if(g<=1||y<=1)return null;const b=c&&r.length>=e*t*3;if(b){if(i<0||i>2)return null}else if(i<0||i>3)return null;let v=0,A=0;if(b)v=p-d,A=_-u;else{for(let f=d;f<p;f++)f%2===(i&1)&&v++;for(let f=u;f<_;f++)f%2===i>>1&&A++}if(v<4||A<4)return null;const T=Math.max(1,h),P=Math.floor(2*v/(T+1)),O=Math.floor(2*A/(T+1));if(P<4||O<4)return null;const E=new He(P),U=new He(O),k=new Float32Array(Math.floor(P/2)+1),W=new Float32Array(Math.floor(O/2)+1);let q=0,X=0;const x=new Float32Array(v),B=new Float32Array(A);for(let f=u;f<_;f++){const C=f%2;if(!b&&i>>1!==C||o.rows&&(b?f&1:f>>1&1)!==o.rowGroup)continue;const D=f*e;let S=0;for(let w=d;w<p;w++){if(!b&&w%2!==(i&1)||S>=v)continue;const R=b?r[(D+w)*3+i]:r[D+w],J=Number.isFinite(R)?R:0;x[S++]=J}for(let w=0;w<T;w++){const R=Math.floor(w*P/2);if(R+P>v)break;E.calculateSpectrumWindow(x,R,k,a),q++}if(!l)break}for(let f=d;f<p;f++){const C=f%2;if(!b&&(i&1)!==C||o.cols&&(b?f&1:f>>1&1)!==o.colGroup)continue;let D=0;for(let S=u;S<_;S++){if(!b&&S%2!==i>>1||D>=A)continue;const w=b?r[(S*e+f)*3+i]:r[S*e+f],R=Number.isFinite(w)?w:0;B[D++]=R}for(let S=0;S<T;S++){const w=Math.floor(S*O/2);if(w+O>A)break;U.calculateSpectrumWindow(B,w,W,a),X++}if(!l)break}if(q>0)for(let f=0;f<k.length;f++)k[f]/=q;if(X>0)for(let f=0;f<W.length;f++)W[f]/=X;return{avgHoriz:Array.from(k),avgVert:Array.from(W),sampleCount:{h:q,v:X}}}const en=192*1024*1024,Ye=new Map;let Ge=0,tn=Promise.resolve();const br=r=>{r.lastAccess=Date.now()},xr=(r=new Set)=>{if(Ge<=en)return;const e=Array.from(Ye.entries()).filter(([t])=>!r.has(t)).sort((t,n)=>t[1].lastAccess-n[1].lastAccess);for(const[t,n]of e){if(Ge<=en)break;Ge-=n.sizeBytes,Ye.delete(t)}},Sr=r=>{var e;return r.data.byteLength+(((e=r.floatData)==null?void 0:e.byteLength)||0)},xn=async r=>{const e=Ye.get(r.key);if(e)return br(e),e.raw;if(!r.file)throw new Error(`Worker cache miss for ${r.key}`);const t=await r.file.arrayBuffer(),n=await yr(t),i={raw:n,sizeBytes:Sr(n),lastAccess:Date.now()};return Ye.set(r.key,i),Ge+=i.sizeBytes,xr(new Set([r.key])),n},kr=async r=>{for(const e of r.files)await xn(e);self.postMessage({id:r.id,type:"result",success:!0})},Cr=async r=>{const e=new Map;for(const i of r.files){const s=await xn(i);for(const a of r.channelIndices){const o=wr(s.data,s.width,s.height,r.roi,a,!1,r.useHanning,r.interlace,r.doAverage,r.nSegments,!!s.isThreePlane);if(!o)continue;let l=e.get(a);if(!l)l={sumH:new Float32Array(o.avgHoriz.length),sumV:new Float32Array(o.avgVert.length),sampleCount:o.sampleCount,contributingFrames:0},e.set(a,l);else if(l.sumH.length!==o.avgHoriz.length||l.sumV.length!==o.avgVert.length)throw new Error("Stack FFT requires matching frame dimensions and channel layout.");for(let h=0;h<o.avgHoriz.length;h++)l.sumH[h]+=o.avgHoriz[h];for(let h=0;h<o.avgVert.length;h++)l.sumV[h]+=o.avgVert[h];l.sampleCount=o.sampleCount,l.contributingFrames+=1}self.postMessage({id:r.id,type:"progress",completed:1})}const t=Array.from(e.entries()).map(([i,s])=>({ch:i,sumH:s.sumH,sumV:s.sumV,sampleCount:s.sampleCount,contributingFrames:s.contributingFrames})),n=[];for(const i of t)n.push(i.sumH.buffer,i.sumV.buffer);self.postMessage({id:r.id,type:"result",success:!0,processedFrames:r.files.length,results:t},n)};self.onmessage=r=>{const e=r.data;tn=tn.then(async()=>{try{e.type==="preload"?await kr(e):await Cr(e)}catch(t){self.postMessage({id:e.id,type:"result",success:!1,error:(t==null?void 0:t.message)||"Power spectrum worker failed"})}})};
