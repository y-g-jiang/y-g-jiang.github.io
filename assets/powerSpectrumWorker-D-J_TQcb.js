var xr=Object.defineProperty;var Sr=(n,e,t)=>e in n?xr(n,e,{enumerable:!0,configurable:!0,writable:!0,value:t}):n[e]=t;var z=(n,e,t)=>Sr(n,typeof e!="symbol"?e+"":e,t);class He{constructor(e){z(this,"size");z(this,"isPowerOfTwo");z(this,"_real");z(this,"_imag");z(this,"_scratch",null);z(this,"_rev",null);z(this,"_m",0);z(this,"_internalFFT",null);z(this,"_chirpReal",null);z(this,"_chirpImag",null);z(this,"_bReal",null);z(this,"_bImag",null);z(this,"_hanning",null);z(this,"_windowSumSq",0);this.size=e,this.isPowerOfTwo=(e&e-1)===0&&e>0,this._real=new Float32Array(e),this._imag=new Float32Array(e),this.isPowerOfTwo?this.initRadix2():this.initBluestein()}initRadix2(){const e=this.size,t=Math.log2(e);this._rev=new Uint32Array(e);for(let r=0;r<e;r++){let i=0,s=r;for(let a=0;a<t;a++)i=i<<1|s&1,s>>>=1;this._rev[r]=i}}initBluestein(){const e=this.size;this._m=Math.pow(2,Math.ceil(Math.log2(2*e-1))),this._internalFFT=new He(this._m),this._chirpReal=new Float32Array(e),this._chirpImag=new Float32Array(e);for(let i=0;i<e;i++){const s=-Math.PI*(i*i)/e;this._chirpReal[i]=Math.cos(s),this._chirpImag[i]=Math.sin(s)}const t=new Float32Array(this._m),r=new Float32Array(this._m);for(let i=0;i<e;i++)t[i]=this._chirpReal[i],r[i]=-this._chirpImag[i];for(let i=1;i<e;i++)t[this._m-i]=t[i],r[this._m-i]=r[i];this._internalFFT.transform(t,r),this._bReal=new Float32Array(this._internalFFT._real),this._bImag=new Float32Array(this._internalFFT._imag)}initHanning(){if(this._hanning)return;const e=this.size;this._hanning=new Float32Array(e);let t=0;for(let r=0;r<e;r++){const i=.5*(1-Math.cos(2*Math.PI*r/(e-1)));this._hanning[r]=i,t+=i*i}this._windowSumSq=t}transform(e,t){this.isPowerOfTwo?this.transformRadix2(e,t):this.transformBluestein(e,t)}transformRadix2(e,t){const r=this.size,i=this._rev,s=this._real,a=this._imag;if(e===s)for(let o=0;o<r;o++){const l=i[o];if(o<l){const h=s[o],c=a[o];s[o]=s[l],a[o]=a[l],s[l]=h,a[l]=c}}else for(let o=0;o<r;o++){const l=i[o];s[o]=e[l],a[o]=t?t[l]:0}for(let o=2;o<=r;o*=2){const l=o/2,h=-2*Math.PI/o,c=Math.cos(h),d=Math.sin(h);for(let u=0;u<r;u+=o){let p=1,_=0;for(let g=0;g<l;g++){const y=u+g,b=u+g+l,v=p*s[b]-_*a[b],A=p*a[b]+_*s[b],T=s[y],P=a[y];s[y]=T+v,a[y]=P+A,s[b]=T-v,a[b]=P-A;const O=p*c-_*d,E=p*d+_*c;p=O,_=E}}}}transformBluestein(e,t){const r=this.size,i=this._m,s=this._internalFFT,a=s._real,o=s._imag;a.fill(0),o.fill(0);for(let d=0;d<r;d++){const u=e[d],p=t?t[d]:0,_=this._chirpReal[d],g=this._chirpImag[d];a[d]=u*_-p*g,o[d]=u*g+p*_}s.transformRadix2(a,o);for(let d=0;d<i;d++){const u=s._real[d],p=s._imag[d],_=this._bReal[d],g=this._bImag[d];s._real[d]=u*_-p*g,s._imag[d]=u*g+p*_}const l=s._real,h=s._imag;for(let d=0;d<i;d++)h[d]=-h[d];s.transformRadix2(l,h);const c=1/i;for(let d=0;d<r;d++){const u=s._real[d]*c,p=-s._imag[d]*c,_=this._chirpReal[d],g=this._chirpImag[d];this._real[d]=u*_-p*g,this._imag[d]=u*g+p*_}}calculateSpectrum(e,t,r=!1){const i=this.size;let s=0;for(let c=0;c<i;c++)s+=e[c];const a=s/i;this._scratch||(this._scratch=new Float32Array(i));const o=this._scratch;if(r){this.initHanning();const c=this._hanning;for(let d=0;d<i;d++)o[d]=(e[d]-a)*c[d]}else for(let c=0;c<i;c++)o[c]=e[c]-a;this.transform(o);const l=t.length;let h=1/i;r&&this._windowSumSq>0&&(h=1/this._windowSumSq);for(let c=0;c<l;c++){const d=this._real[c],u=this._imag[c];t[c]+=(d*d+u*u)*h}}calculateSpectrumWindow(e,t,r,i=!1){const s=this.size;let a=0;for(let d=0;d<s;d++)a+=e[t+d];const o=a/s;this._scratch||(this._scratch=new Float32Array(s));const l=this._scratch;if(i){this.initHanning();const d=this._hanning;for(let u=0;u<s;u++)l[u]=(e[t+u]-o)*d[u]}else for(let d=0;d<s;d++)l[d]=e[t+d]-o;this.transform(l);const h=r.length;let c=1/s;i&&this._windowSumSq>0&&(c=1/this._windowSumSq);for(let d=0;d<h;d++){const u=this._real[d],p=this._imag[d];r[d]+=(u*u+p*p)*c}}}const kr={"Sony ILCE-7RM5":"0.82 -0.2976 -0.0719 -0.4296 1.2053 0.2532 -0.0429 0.1282 0.5774"};let at=null;async function Cr(n){return at||(at=(async()=>{if(typeof window.loadPyodide!="function")throw new Error("Pyodide missing: window.loadPyodide not found.");const e=await window.loadPyodide();return await e.loadPackage("numpy"),e})()),at}let ot=null,Rt=!1;async function bt(){var t;ot||(ot=(async()=>{const r=await import("./joraw2-Bb3_vNP4.js");if(typeof r.default!="function")throw new Error("JoRaw2 WASM import failed");const i=new URL("/assets/joraw2-3YkywkGx.wasm",import.meta.url).href;return r.default({locateFile(s,a){return s.endsWith("joraw2.wasm")?i:a+s}})})());const e=(await ot).LibRaw;if(!e)throw new Error("JoRaw2 class not found");if(!Rt){const r=new e;try{if(typeof r.runtimeInfo!="function")throw new Error("JoRaw2 runtime identity is missing");const i=r.runtimeInfo();if((i==null?void 0:i.wrapper)!=="joraw2"||!String((i==null?void 0:i.librawVersion)||"").startsWith("0.22.2")||!(i!=null&&i.nikonHe)||!(i!=null&&i.nikonHeStar))throw new Error(`Unexpected JoRaw2 runtime: ${JSON.stringify(i)}`);Rt=!0}finally{(t=r.delete)==null||t.call(r)}}return e}var vr=`#!/usr/bin/env python3\r
"""Container parser for the Sony ARW6/LLVC3 raw strip.\r
\r
This stops at packet records and control fields. Coefficient entropy decoding\r
lives in llvc3_entropy.py; this file is mostly the thing I dump to JSON and diff\r
against Imaging Edge traces when the packet framing looks suspicious.\r
"""\r
\r
from __future__ import annotations\r
\r
import argparse\r
import json\r
import struct\r
from dataclasses import asdict, dataclass\r
from pathlib import Path\r
from typing import Any\r
\r
\r
RAW_STREAM_OFFSET = 0x200\r
\r
\r
def u16be(buf: bytes, off: int) -> int:\r
    return struct.unpack_from(">H", buf, off)[0]\r
\r
\r
def u16le(buf: bytes, off: int) -> int:\r
    return struct.unpack_from("<H", buf, off)[0]\r
\r
\r
def u32le(buf: bytes, off: int) -> int:\r
    return struct.unpack_from("<I", buf, off)[0]\r
\r
\r
def u32be(buf: bytes, off: int) -> int:\r
    return struct.unpack_from(">I", buf, off)[0]\r
\r
\r
class BitReader:\r
    """Small MSB-first reader for packet headers."""\r
\r
    def __init__(self, data: bytes) -> None:\r
        self.data = data\r
        self.pos = 0\r
\r
    def read(self, nbits: int) -> int:\r
        if nbits < 0:\r
            raise ValueError("negative bit count")\r
        out = 0\r
        for _ in range(nbits):\r
            if self.pos >= len(self.data) * 8:\r
                raise EOFError("packet bitstream exhausted")\r
            byte = self.data[self.pos >> 3]\r
            bit = (byte >> (7 - (self.pos & 7))) & 1\r
            out = (out << 1) | bit\r
            self.pos += 1\r
        return out\r
\r
\r
@dataclass\r
class TiffRawInfo:\r
    ifd_offset: int\r
    width: int\r
    height: int\r
    bits_per_sample: int\r
    compression: int\r
    photometric: int\r
    strip_offset: int\r
    strip_byte_count: int\r
    cfa_pattern: list[int] | None\r
    black_level_tag_0x7310: list[int] | None\r
    white_level: int | None\r
    default_crop_origin: list[int] | None\r
    default_crop_size: list[int] | None\r
\r
\r
@dataclass\r
class LlvcHeader:\r
    magic: str\r
    sequence_or_version: int\r
    coded_width: int\r
    coded_half_height: int\r
    logical_height: int\r
    decoded_bits: int\r
    component_count: int\r
    mode: int\r
    flags_low10: int\r
\r
\r
@dataclass\r
class LlvcStreamInfo:\r
    index: int\r
    offset: int\r
    length: int\r
    header: LlvcHeader\r
    tile_x: int\r
    tile_y: int\r
    tile_width: int\r
    tile_height: int\r
\r
\r
@dataclass\r
class DirectoryEntry:\r
    group: int\r
    index: int\r
    start: int\r
    length: int\r
    mode_in_decoder_trace: int | None = None\r
\r
\r
@dataclass\r
class PacketRecord:\r
    index: int\r
    byte_length: int\r
    selectors: list[int]\r
    payload_offset: int\r
\r
\r
@dataclass\r
class PacketHeader:\r
    group: int\r
    index: int\r
    stream_offset: int\r
    directory_length: int\r
    control_count: int\r
    extra_count: int\r
    tag4: int\r
    reserved4: int\r
    type2: int\r
    control_words: int\r
    block_count: int\r
    width_marker: int\r
    reserved8: int\r
    skipped_u8: list[int]\r
    header_bits: int\r
    control_bytes: int\r
    total_bytes: int\r
    payload_bytes_from_records: int\r
    validation: dict[str, bool]\r
    records: list[PacketRecord]\r
    first_records: list[PacketRecord]\r
    last_records: list[PacketRecord]\r
\r
\r
def parse_tiff_value(data: bytes, bo: str, typ: int, cnt: int, raw: bytes) -> Any:\r
    fmt_by_type = {1: "B", 3: "H", 4: "I", 8: "h", 9: "i"}\r
    if typ in fmt_by_type:\r
        fmt = bo + fmt_by_type[typ] * cnt\r
        size = struct.calcsize(fmt)\r
        val = list(struct.unpack(fmt, raw[:size]))\r
        return val[0] if len(val) == 1 else val\r
    if typ in (2, 7):\r
        return list(raw[:cnt])\r
    return raw.hex(" ")\r
\r
\r
def iter_tiff_ifds(data: bytes) -> list[tuple[int, dict[int, Any]]]:\r
    if data[:2] == b"II":\r
        bo = "<"\r
    elif data[:2] == b"MM":\r
        bo = ">"\r
    else:\r
        raise ValueError("input is not a TIFF/ARW file")\r
\r
    def u16(off: int) -> int:\r
        return struct.unpack_from(bo + "H", data, off)[0]\r
\r
    def u32(off: int) -> int:\r
        return struct.unpack_from(bo + "I", data, off)[0]\r
\r
    type_sizes = {1: 1, 2: 1, 3: 2, 4: 4, 5: 8, 7: 1, 8: 2, 9: 4, 10: 8}\r
    stack = [u32(4)]\r
    seen: set[int] = set()\r
    out: list[tuple[int, dict[int, Any]]] = []\r
    while stack:\r
        off = stack.pop()\r
        if off in seen or off <= 0 or off >= len(data):\r
            continue\r
        seen.add(off)\r
        n = u16(off)\r
        tags: dict[int, Any] = {}\r
        for i in range(n):\r
            ent = off + 2 + i * 12\r
            tag = u16(ent)\r
            typ = u16(ent + 2)\r
            cnt = u32(ent + 4)\r
            value_area = ent + 8\r
            size = type_sizes.get(typ, 1) * cnt\r
            raw = data[value_area : value_area + 4] if size <= 4 else data[u32(value_area) : u32(value_area) + size]\r
            tags[tag] = parse_tiff_value(data, bo, typ, cnt, raw)\r
            if tag == 0x014A:\r
                offsets = tags[tag] if isinstance(tags[tag], list) else [tags[tag]]\r
                stack.extend(int(x) for x in offsets)\r
        next_ifd = u32(off + 2 + n * 12)\r
        if next_ifd:\r
            stack.append(next_ifd)\r
        out.append((off, tags))\r
    return out\r
\r
\r
def scalar_or_first(value: Any) -> int:\r
    if isinstance(value, list):\r
        return int(value[0])\r
    return int(value)\r
\r
\r
def find_raw_subifd(path: Path) -> tuple[TiffRawInfo, bytes]:\r
    data = path.read_bytes()\r
    for off, tags in iter_tiff_ifds(data):\r
        if tags.get(0x0103) == 32766 and tags.get(0x0106) == 32803:\r
            strip_offset = scalar_or_first(tags[0x0111])\r
            strip_len = scalar_or_first(tags[0x0117])\r
            info = TiffRawInfo(\r
                ifd_offset=off,\r
                width=int(tags[0x0100]),\r
                height=int(tags[0x0101]),\r
                bits_per_sample=int(tags[0x0102]),\r
                compression=int(tags[0x0103]),\r
                photometric=int(tags[0x0106]),\r
                strip_offset=strip_offset,\r
                strip_byte_count=strip_len,\r
                cfa_pattern=list(tags.get(0x828E, [])) or None,\r
                black_level_tag_0x7310=list(tags.get(0x7310, [])) if isinstance(tags.get(0x7310), list) else None,\r
                white_level=int(tags[0xC61D]) if 0xC61D in tags else None,\r
                default_crop_origin=list(tags.get(0xC61F, [])) if isinstance(tags.get(0xC61F), list) else None,\r
                default_crop_size=list(tags.get(0xC620, [])) if isinstance(tags.get(0xC620), list) else None,\r
            )\r
            return info, data[strip_offset : strip_offset + strip_len]\r
    raise ValueError("no ARW6 LLVC raw SubIFD found")\r
\r
\r
def parse_llvc_header(stream: bytes) -> LlvcHeader:\r
    word_c = u16be(stream, 0x0C)\r
    word_e = u16be(stream, 0x0E)\r
    return LlvcHeader(\r
        magic=stream[:4].decode("ascii", "replace"),\r
        sequence_or_version=u32le(stream, 0x04),\r
        coded_width=u16be(stream, 0x08),\r
        coded_half_height=u16be(stream, 0x0A),\r
        logical_height=u16be(stream, 0x0A) * 2,\r
        decoded_bits=(word_c >> 4) & 0x3F,\r
        component_count=word_e >> 13,\r
        mode=(word_e >> 10) & 0x03,\r
        flags_low10=word_e & 0x03FF,\r
    )\r
\r
\r
def initial_group_lengths(stream: bytes) -> list[int]:\r
    return [((u32be(stream, 0x10 + off) >> 4) & 0x0FFFFFF0) for off in (0, 3, 6, 9, 12)]\r
\r
\r
def llvc_stream_length(stream: bytes) -> int:\r
    """Return the byte span occupied by one LLVC3 stream."""\r
\r
    return 0x80 + sum(initial_group_lengths(stream))\r
\r
\r
def is_plausible_llvc_header(header: LlvcHeader) -> bool:\r
    return (\r
        header.magic in {"A000", "0000"}\r
        and header.coded_width > 0\r
        and header.coded_half_height > 0\r
        and header.decoded_bits == 16\r
        and header.component_count == 3\r
        and header.mode == 3\r
    )\r
\r
\r
def find_llvc_streams(strip: bytes) -> list[LlvcStreamInfo]:\r
    """Find all LLVC3 streams inside an ARW6 raw strip."""\r
\r
    streams: list[LlvcStreamInfo] = []\r
    count = u32le(strip, 0) if len(strip) >= 4 else 0\r
    if 1 <= count <= 16 and len(strip) >= RAW_STREAM_OFFSET + 0x80:\r
        for index in range(count):\r
            entry = 0x08 + index * 0x18\r
            if entry + 0x18 > len(strip):\r
                streams = []\r
                break\r
            table_offset = u32le(strip, entry)\r
            tile_x = u32le(strip, entry + 0x08)\r
            tile_y = u32le(strip, entry + 0x0C)\r
            tile_width = u32le(strip, entry + 0x10)\r
            tile_height = u32le(strip, entry + 0x14)\r
            pos = table_offset if table_offset else RAW_STREAM_OFFSET\r
            found_pos: int | None = None\r
            search_end = min(len(strip) - 0x80, pos + 0x1000)\r
            for candidate in range(pos, search_end + 1, 0x10):\r
                try:\r
                    candidate_header = parse_llvc_header(strip[candidate:])\r
                except Exception:\r
                    continue\r
                if is_plausible_llvc_header(candidate_header):\r
                    found_pos = candidate\r
                    header = candidate_header\r
                    break\r
            if found_pos is None:\r
                streams = []\r
                break\r
            pos = found_pos\r
            if not is_plausible_llvc_header(header):\r
                streams = []\r
                break\r
            length = llvc_stream_length(strip[pos:])\r
            if length <= 0x80 or pos + length > len(strip):\r
                streams = []\r
                break\r
            streams.append(\r
                LlvcStreamInfo(\r
                    index=index,\r
                    offset=pos,\r
                    length=length,\r
                    header=header,\r
                    tile_x=tile_x,\r
                    tile_y=tile_y,\r
                    tile_width=tile_width or header.coded_width,\r
                    tile_height=tile_height or header.logical_height,\r
                )\r
            )\r
        if len(streams) == count:\r
            return streams\r
\r
    for pos in range(RAW_STREAM_OFFSET, max(RAW_STREAM_OFFSET, len(strip) - 0x80), 0x10):\r
        try:\r
            header = parse_llvc_header(strip[pos:])\r
        except Exception:\r
            continue\r
        if not is_plausible_llvc_header(header):\r
            continue\r
        try:\r
            length = llvc_stream_length(strip[pos:])\r
        except Exception:\r
            continue\r
        if length <= 0x80 or pos + length > len(strip):\r
            continue\r
        streams.append(\r
            LlvcStreamInfo(\r
                index=len(streams),\r
                offset=pos,\r
                length=length,\r
                header=header,\r
                tile_x=0,\r
                tile_y=0,\r
                tile_width=header.coded_width,\r
                tile_height=header.logical_height,\r
            )\r
        )\r
    return streams\r
\r
\r
def parse_directory(stream: bytes) -> tuple[int, list[DirectoryEntry], list[dict[str, Any]]]:\r
    group_lengths = initial_group_lengths(stream)\r
    pos = 0x30\r
    base = 0\r
    entries: list[DirectoryEntry] = []\r
    groups: list[dict[str, Any]] = []\r
    for group, group_len in enumerate(group_lengths):\r
        n_entries = stream[pos] & 0x0F\r
        values = [(u32be(stream, pos + off) & 0x00FFFFFF) << 4 for off in (0, 3, 6, 9, 12)]\r
        consumed = 0x10\r
        if n_entries >= 5:\r
            values.extend((u32be(stream, pos + 0x10 + off) & 0x00FFFFFF) << 4 for off in (0, 3, 6, 9))\r
            consumed = 0x20\r
        local = 0\r
        for index in range(n_entries):\r
            entries.append(DirectoryEntry(group=group, index=index, start=0x80 + base + local, length=values[index]))\r
            local += values[index]\r
        groups.append(\r
            {\r
                "group": group,\r
                "declared_length": group_len,\r
                "entry_count": n_entries,\r
                "entry_lengths": values[:n_entries],\r
                "entry_sum": local,\r
                "directory_offset": pos,\r
                "directory_bytes": consumed,\r
                "sum_matches_declared": local == group_len,\r
            }\r
        )\r
        base += group_len\r
        pos += consumed\r
    return pos, entries, groups\r
\r
\r
def parse_packet(stream: bytes, entry: DirectoryEntry) -> PacketHeader:\r
    packet = stream[entry.start : entry.start + entry.length]\r
    br = BitReader(packet)\r
    control_count = br.read(16)\r
    extra_count = br.read(24)\r
    tag4 = br.read(4)\r
    reserved4 = br.read(4)\r
    type2 = br.read(2)\r
    control_words = br.read(6)\r
    block_count = br.read(16)\r
    width_marker = br.read(8)\r
    reserved8 = br.read(8)\r
    skipped = [br.read(8) for _ in range(5)]\r
    control_bytes = (control_count + 1) << 4\r
    total_bytes = (control_count + 1 + extra_count) << 4\r
    records: list[PacketRecord] = []\r
    cursor = control_bytes\r
    for i in range(block_count):\r
        byte_len = br.read(16)\r
        selectors = [br.read(4) for _ in range(type2)]\r
        records.append(PacketRecord(i, byte_len, selectors, cursor))\r
        cursor += byte_len\r
    # Native code consumes the 6-bit field while walking the header, then the\r
    # validation at 0x1aa717 compares the original 16-bit control_count with:\r
    # ceil((type2 + 4) * block_count * 4 / 128).\r
    formula = ((type2 + 4) * block_count * 4 + 0x7F) >> 7\r
    validation = {\r
        "directory_length_matches_total": entry.length == total_bytes,\r
        "tag4_is_4": tag4 == 4,\r
        "reserved4_is_0": reserved4 == 0,\r
        "type2_is_1_or_3": type2 in (1, 3),\r
        "reserved6_is_0": control_words == 0,\r
        "control_count_formula": control_count == formula,\r
        "block_count_le_300": block_count <= 300,\r
        "width_marker_is_0x10": width_marker == 0x10,\r
        "reserved8_is_0": reserved8 == 0,\r
        "skipped_bytes_are_0": all(x == 0 for x in skipped),\r
        "record_payload_within_total": cursor <= total_bytes,\r
    }\r
    return PacketHeader(\r
        group=entry.group,\r
        index=entry.index,\r
        stream_offset=entry.start,\r
        directory_length=entry.length,\r
        control_count=control_count,\r
        extra_count=extra_count,\r
        tag4=tag4,\r
        reserved4=reserved4,\r
        type2=type2,\r
        control_words=control_words,\r
        block_count=block_count,\r
        width_marker=width_marker,\r
        reserved8=reserved8,\r
        skipped_u8=skipped,\r
        header_bits=br.pos,\r
        control_bytes=control_bytes,\r
        total_bytes=total_bytes,\r
        payload_bytes_from_records=cursor,\r
        validation=validation,\r
        records=records,\r
        first_records=records[:8],\r
        last_records=records[-4:],\r
    )\r
\r
\r
def derive_metrics(raw_info: TiffRawInfo, strip_len: int, stream_len: int) -> dict[str, Any]:\r
    samples = raw_info.width * raw_info.height\r
    return {\r
        "samples": samples,\r
        "decoded_u16_bytes": samples * 2,\r
        "decoded_14bit_packed_bytes": samples * 14 / 8,\r
        "strip_bytes_per_sample": strip_len / samples,\r
        "strip_bits_per_sample": strip_len * 8 / samples,\r
        "stream_bits_per_sample_excluding_0x200_preamble": stream_len * 8 / samples,\r
        "compression_ratio_vs_u16": (samples * 2) / strip_len,\r
        "compression_ratio_vs_14bit_packed": (samples * 14 / 8) / strip_len,\r
    }\r
\r
\r
def main() -> None:\r
    ap = argparse.ArgumentParser()\r
    ap.add_argument("path", nargs="?", default="DSC00089.ARW")\r
    ap.add_argument("--out", default="out/reverse/llvc3_bitstream_probe.json")\r
    ns = ap.parse_args()\r
\r
    raw_info, strip = find_raw_subifd(Path(ns.path))\r
    stream = strip[RAW_STREAM_OFFSET:]\r
    header = parse_llvc_header(stream)\r
    consumed, entries, groups = parse_directory(stream)\r
    packets = [parse_packet(stream, entry) for entry in entries]\r
    all_valid = all(all(p.validation.values()) for p in packets)\r
    type_counts: dict[str, int] = {}\r
    for p in packets:\r
        type_counts[str(p.type2)] = type_counts.get(str(p.type2), 0) + 1\r
    result = {\r
        "input": str(ns.path),\r
        "raw_subifd": asdict(raw_info),\r
        "llvc_header": asdict(header),\r
        "stream_offset_inside_raw_strip": RAW_STREAM_OFFSET,\r
        "strip_preamble_first_32": strip[:32].hex(" "),\r
        "directory": {\r
            "consumed_bytes_after_stream_header": consumed - 0x10,\r
            "packet_base_offset": 0x80,\r
            "groups": groups,\r
            "entries": [asdict(e) for e in entries],\r
        },\r
        "packets": [asdict(p) for p in packets],\r
        "summary": {\r
            "packet_count": len(packets),\r
            "packet_type_counts": type_counts,\r
            "all_packet_validations_pass": all_valid,\r
            "total_packet_bytes": sum(p.directory_length for p in packets),\r
            "payload_bytes_from_records": sum(p.payload_bytes_from_records for p in packets),\r
            "metrics": derive_metrics(raw_info, len(strip), len(stream)),\r
        },\r
    }\r
    out = Path(ns.out)\r
    out.parent.mkdir(parents=True, exist_ok=True)\r
    out.write_text(json.dumps(result, indent=2), encoding="utf-8")\r
    print(json.dumps(result["summary"], indent=2))\r
\r
\r
if __name__ == "__main__":\r
    main()\r
`,Pr=`#!/usr/bin/env python3\r
"""Entropy-side notes for the Sony ARW6/LLVC3 stream.\r
\r
The bit reader and 4-lane coefficient paths here came straight out of Imaging\r
Edge traces. I kept the code narrow on purpose: first replay a row, then a\r
packet, then let the higher-level decoder stitch the pieces together.\r
"""\r
\r
from __future__ import annotations\r
\r
import argparse\r
import json\r
from dataclasses import dataclass\r
from pathlib import Path\r
from typing import Iterable\r
\r
from llvc3_bitstream_probe import find_llvc_streams, find_raw_subifd\r
\r
\r
@dataclass\r
class NativeBitReader:\r
    """MSB-first 64-bit reader with the same state layout Sony uses.\r
\r
    The offsets are worth keeping close by:\r
\r
      +0x00 current big-endian 64-bit word\r
      +0x08 pointer to the last loaded word\r
      +0x10 bit position / remaining bits in current word\r
      +0x14 64-bit words left to load\r
      +0x1c status\r
\r
    The pointer is initialized to eight bytes before the record payload. On\r
    underflow the decoder advances by eight bytes and loads a big-endian word.\r
    """\r
\r
    data: bytes\r
    ptr: int = -8\r
    cur: int = 0\r
    bit: int = 0\r
    words_left: int = 0\r
    status: int = 0\r
\r
    @classmethod\r
    def for_record(cls, packet: bytes, payload_offset: int, byte_length: int) -> "NativeBitReader":\r
        words = (((byte_length + 7) // 8) + 1) & ~1\r
        return cls(data=packet, ptr=payload_offset - 8, words_left=words)\r
\r
    def _load_next_word(self) -> int:\r
        self.words_left -= 1\r
        if self.words_left < 0:\r
            self.status = 2\r
            return 0\r
        self.ptr += 8\r
        chunk = self.data[self.ptr : self.ptr + 8]\r
        if len(chunk) < 8:\r
            chunk = chunk + b"\\x00" * (8 - len(chunk))\r
        self.cur = int.from_bytes(chunk, "big")\r
        self.bit = 64\r
        return self.cur\r
\r
    def read_bits(self, nbits: int) -> int:\r
        if nbits == 0:\r
            return 0\r
        if nbits < 0:\r
            raise ValueError("negative bit count")\r
        out = 0\r
        remaining = nbits\r
        while remaining > 0:\r
            if self.bit <= 0:\r
                self._load_next_word()\r
                if self.status:\r
                    return out << remaining\r
            take = min(remaining, self.bit)\r
            self.bit -= take\r
            out = (out << take) | ((self.cur >> self.bit) & ((1 << take) - 1))\r
            remaining -= take\r
        return out\r
\r
    def skip_zero_words_to_one(self, initial_bit_count: int) -> int:\r
        """Count leading zeros until the next one, in the 0x1a9080 style."""\r
\r
        total = 0\r
        while True:\r
            if self.bit <= 0:\r
                self._load_next_word()\r
                if self.status:\r
                    return total\r
            if self.cur == 0:\r
                total += self.bit\r
                self.bit = 0\r
                continue\r
            # Count zeros among currently available MSB-side bits.\r
            while self.bit > 0:\r
                self.bit -= 1\r
                if (self.cur >> self.bit) & 1:\r
                    return total\r
                total += 1\r
                if total >= initial_bit_count and initial_bit_count > 0:\r
                    return total\r
\r
    def read_unary_zeros_plus_one(self) -> int:\r
        zeros = 0\r
        while not self.status:\r
            bit = self.read_bits(1)\r
            if self.status:\r
                break\r
            if bit:\r
                return zeros + 1\r
            zeros += 1\r
        return zeros + 1\r
\r
\r
def split_lanes(packed: int, width: int) -> list[int]:\r
    mask = (1 << width) - 1\r
    return [\r
        (packed >> (3 * width)) & mask,\r
        (packed >> (2 * width)) & mask,\r
        (packed >> width) & mask,\r
        packed & mask,\r
    ]\r
\r
\r
def magnitude4(br: NativeBitReader, width: int, shift: int) -> list[int]:\r
    """Magnitude expansion from 0x1a8b00."""\r
\r
    if width <= 0:\r
        raw = [0, 0, 0, 0]\r
    else:\r
        raw = split_lanes(br.read_bits(width * 4), width)\r
    if shift <= 0:\r
        return raw\r
    out: list[int] = []\r
    for x in raw:\r
        if x <= 0:\r
            out.append(0)\r
        else:\r
            # Visible SIMD shape: ((2*x + 1) << (shift - 1)) - (x & 1).\r
            out.append(((2 * x + 1) << (shift - 1)) - (x & 1))\r
    return out\r
\r
\r
def apply_sign4(br: NativeBitReader, coeffs: Iterable[int]) -> list[int]:\r
    """Apply one sign bit to each positive lane; see 0x1a8dd0."""\r
\r
    out: list[int] = []\r
    for x in coeffs:\r
        if x > 0:\r
            bit = br.read_bits(1)\r
            out.append(x - 2 * x * bit)\r
        else:\r
            out.append(x)\r
    return out\r
\r
\r
def update_width(br: NativeBitReader, width: int) -> int:\r
    """Adaptive width update shared by 0x1a9080 and 0x1ac060.\r
\r
    Prefix structure from traces and branch shape:\r
\r
      0      keep width\r
      10 U   increase by unary U, where U is zeros+1 terminated by a one\r
      11 U   decrease by unary U, clipped at zero\r
    """\r
\r
    if br.read_bits(1) == 0 or br.status:\r
        return width\r
    if br.read_bits(1) == 0 or br.status:\r
        return width + br.read_unary_zeros_plus_one()\r
\r
    # Decrement path: if the run would cross zero, native code consumes only\r
    # width-1 zeros and leaves the terminating one for the next state; that is\r
    # the 0x1ac4ff/0x1ac50b saturation branch.\r
    for zeros in range(max(0, width - 1)):\r
        if br.read_bits(1):\r
            return width - (zeros + 1)\r
    return 0\r
\r
\r
def read_initial_width(br: NativeBitReader) -> int:\r
    """Read the row's first adaptive width with the 0x1a9080 zero state."""\r
\r
    return update_width(br, 0)\r
\r
\r
def read_zero_run(br: NativeBitReader, remaining: int) -> int:\r
    """Read the zero-group run used when the adaptive width is zero."""\r
\r
    if remaining <= 1:\r
        return remaining\r
    max_prefix = (remaining - 1).bit_length()\r
    zeros = 0\r
    while zeros < max_prefix:\r
        if br.read_bits(1):\r
            break\r
        if br.status:\r
            return remaining\r
        zeros += 1\r
    else:\r
        return remaining\r
    base = 1 << zeros\r
    if base >= remaining:\r
        return remaining\r
    extra = br.read_bits(zeros) if zeros else 0\r
    run = base + extra\r
    return min(run, remaining)\r
\r
\r
def packet_row_multiplier(group: int) -> int:\r
    """Number of output rows produced by one packet record for this scale."""\r
\r
    if group < 0:\r
        raise ValueError("negative packet group")\r
    if group == 4:\r
        return 8\r
    if group == 0:\r
        return 1\r
    return 1 << (group - 1)\r
\r
\r
def read_width_after_zero_run(br: NativeBitReader) -> int:\r
    """Positive width code after a zero run: leading zeros plus the one bit."""\r
\r
    return br.read_unary_zeros_plus_one()\r
\r
\r
def alt4(br: NativeBitReader, width: int, shift: int, has_next: bool = True) -> tuple[list[int], int]:\r
    """Current reconstruction of the alternate 0x1ac060 4-lane path.\r
\r
    It shares the magnitude transform with 0x1a8b00, then signs only positive\r
    lanes. The width predictor agrees with the small-width traces I have, but\r
    wider rows still deserve spot checks.\r
    """\r
\r
    coeffs = magnitude4(br, width, shift)\r
    next_width = width\r
    if has_next:\r
        next_width = update_width(br, width)\r
    return apply_sign4(br, coeffs), next_width\r
\r
\r
def decode_record_component(\r
    br: NativeBitReader, groups: int, shift: int = 0\r
) -> tuple[list[int], list[int], int]:\r
    """Decode one component from the current bitreader position.\r
\r
    \`groups\` is the number of 4-lane coefficient groups. The packet selector\r
    nibble is passed as \`shift\`; the initial adaptive width is read from the\r
    record payload itself by the 0x1a9080 prefix reader.\r
    """\r
\r
    width = read_initial_width(br)\r
    initial_width = width\r
    coeffs: list[int] = []\r
    widths: list[int] = []\r
    gi = 0\r
    while gi < groups:\r
        if br.status or width > 0x13:\r
            br.status = br.status or 1\r
            coeffs.extend([0, 0, 0, 0] * (groups - gi))\r
            widths.extend([width] * (groups - gi))\r
            break\r
        if width == 0:\r
            run = read_zero_run(br, groups - gi)\r
            coeffs.extend([0, 0, 0, 0] * run)\r
            widths.extend([0] * run)\r
            gi += run\r
            if gi >= groups:\r
                break\r
            width = read_width_after_zero_run(br)\r
            continue\r
        vals, width = alt4(br, width, shift=shift, has_next=gi + 1 < groups)\r
        coeffs.extend(vals)\r
        widths.append(width)\r
        gi += 1\r
    return coeffs, widths, initial_width\r
\r
\r
def load_packet(arw: Path, group: int, index: int, stream_index: int = 0) -> tuple[bytes, dict]:\r
    from llvc3_bitstream_probe import parse_directory, parse_packet\r
\r
    raw_info, strip = find_raw_subifd(arw)\r
    streams = find_llvc_streams(strip)\r
    if not streams:\r
        raise ValueError("no LLVC3 stream found in ARW6 raw strip")\r
    if stream_index < 0 or stream_index >= len(streams):\r
        raise ValueError(f"stream_index {stream_index} out of range for {len(streams)} LLVC3 streams")\r
    stream_info = streams[stream_index]\r
    stream = strip[stream_info.offset : stream_info.offset + stream_info.length]\r
    _consumed, entries, _groups = parse_directory(stream)\r
    entry = next(e for e in entries if e.group == group and e.index == index)\r
    packet_info = parse_packet(stream, entry)\r
    packet = stream[entry.start : entry.start + entry.length]\r
    info = json.loads(json.dumps(packet_info, default=lambda o: o.__dict__))\r
    header = stream_info.header\r
    info["raw_width"] = raw_info.width\r
    info["raw_height"] = raw_info.height\r
    info["work_width"] = header.coded_width\r
    info["work_height"] = header.logical_height\r
    info["stream_index"] = stream_index\r
    info["stream_offset"] = stream_info.offset\r
    info["stream_length"] = stream_info.length\r
    info["llvc_header"] = header.__dict__\r
    return packet, info\r
\r
\r
def replay_row(arw: Path, group: int, index: int, row: int, groups: int = 16) -> dict[str, object]:\r
    packet, info = load_packet(arw, group, index)\r
    rec = info["records"][row]\r
    coeffs, widths, br, initial_width = decode_type1_row(packet, rec, groups)\r
    return {\r
        "packet": {"group": group, "index": index},\r
        "row": row,\r
        "record": rec,\r
        "selector": rec["selectors"][0] if rec["selectors"] else 0,\r
        "initial_width": initial_width,\r
        "coeffs": coeffs,\r
        "next_widths": widths,\r
        "bitreader": {"ptr": br.ptr, "bit": br.bit, "words_left": br.words_left, "status": br.status},\r
    }\r
\r
\r
def decode_type1_row(packet: bytes, rec: dict, groups: int) -> tuple[list[int], list[int], NativeBitReader, int]:\r
    """Decode one type-1 record into 4-lane signed coefficients."""\r
\r
    br = NativeBitReader.for_record(packet, rec["payload_offset"], rec["byte_length"])\r
    shift = rec["selectors"][0] if rec["selectors"] else 0\r
    coeffs, widths, initial_width = decode_record_component(br, groups, shift)\r
    return coeffs, widths, br, initial_width\r
\r
\r
def decode_record_components(\r
    packet: bytes, rec: dict, groups: int, components: int, row_multiplier: int = 1\r
) -> tuple[list[list[int]], list[dict[str, int]]]:\r
    """Decode all components stored in one packet record.\r
\r
    Type-3 records share one payload bitreader across three component streams.\r
    The native block decoder calls the same entropy routine three times with\r
    selector nibbles from the control record, so the reader stays live between\r
    components.\r
    """\r
\r
    br = NativeBitReader.for_record(packet, rec["payload_offset"], rec["byte_length"])\r
    rows: list[list[int]] = []\r
    states: list[dict[str, int]] = []\r
    selectors = rec["selectors"] or []\r
    for ci in range(components):\r
        shift = selectors[ci] if ci < len(selectors) else 0\r
        for _ in range(row_multiplier):\r
            if rec["byte_length"] <= 0:\r
                coeffs = [0] * (groups * 4)\r
            else:\r
                coeffs, _widths, _initial = decode_record_component(br, groups, shift)\r
            rows.append(coeffs)\r
            states.append({"ptr": br.ptr, "bit": br.bit, "words_left": br.words_left, "status": br.status})\r
    return rows, states\r
\r
\r
def infer_packet_width(group: int, packet_type: int, raw_width: int = 7040) -> int:\r
    """Entropy row width in coefficients for this ARW6 sample's packet groups.\r
\r
    Type-3 packets carry wavelet detail subbands. Their entropy rows are\r
    half-width relative to the reconstructed scale reported by the block\r
    object; the horizontal synthesis stage doubles that later.\r
    """\r
\r
    if packet_type == 1:\r
        return raw_width // 2 if group == 4 else raw_width // 16\r
    if packet_type == 3:\r
        if 1 <= group <= 3:\r
            return raw_width // (1 << (5 - group))\r
    raise ValueError(f"width inference not yet known for group {group}, type {packet_type}")\r
\r
\r
def decode_packet_components(\r
    arw: Path, group: int, index: int, out_prefix: Path | None = None, stream_index: int = 0\r
) -> dict[str, object]:\r
    """Decode a type-1 or type-3 packet into one or three int32 component arrays."""\r
\r
    packet, info = load_packet(arw, group, index, stream_index=stream_index)\r
    packet_type = info["type2"]\r
    components = 1 if packet_type == 1 else 3\r
    width = infer_packet_width(group, packet_type, int(info.get("work_width", info.get("raw_width", 7040))))\r
    row_multiplier = packet_row_multiplier(group)\r
    groups_per_row = (width + 3) // 4\r
    planes: list[list[list[int]]] = [[] for _ in range(components)]\r
    final_states: list[list[dict[str, int]]] = []\r
    for rec in info["records"]:\r
        comp_rows, states = decode_record_components(packet, rec, groups_per_row, components, row_multiplier)\r
        for ci in range(components):\r
            start = ci * row_multiplier\r
            end = start + row_multiplier\r
            planes[ci].extend(row[:width] for row in comp_rows[start:end])\r
        final_states.append(states)\r
    outs: list[str] = []\r
    if out_prefix:\r
        import numpy as np\r
\r
        out_prefix.parent.mkdir(parents=True, exist_ok=True)\r
        for ci, rows in enumerate(planes):\r
            arr = np.asarray(rows, dtype=np.int32)\r
            path = out_prefix.with_name(f"{out_prefix.name}_c{ci}.bin")\r
            arr.tofile(path)\r
            outs.append(str(path))\r
    return {\r
        "packet": {"group": group, "index": index, "type": packet_type},\r
        "shape": [len(planes[0]), width],\r
        "row_multiplier": row_multiplier,\r
        "components": components,\r
        "first_row_first16": [plane[0][:16] for plane in planes],\r
        "last_nonempty_row_first16": [\r
            next((plane[i][:16] for i in range(len(plane) - 1, -1, -1) if any(plane[i])), []) for plane in planes\r
        ],\r
        "final_states_tail": final_states[-5:],\r
        "outs": outs,\r
    }\r
\r
\r
def decode_packet_arrays(arw: Path, group: int, index: int, stream_index: int = 0) -> tuple[list["object"], dict[str, object]]:\r
    """Decode a packet and return its int32 component arrays in memory."""\r
\r
    import numpy as np\r
\r
    packet, info = load_packet(arw, group, index, stream_index=stream_index)\r
    packet_type = info["type2"]\r
    components = 1 if packet_type == 1 else 3\r
    width = infer_packet_width(group, packet_type, int(info.get("work_width", info.get("raw_width", 7040))))\r
    row_multiplier = packet_row_multiplier(group)\r
    groups_per_row = (width + 3) // 4\r
    planes: list[list[list[int]]] = [[] for _ in range(components)]\r
    final_states: list[list[dict[str, int]]] = []\r
    for rec in info["records"]:\r
        comp_rows, states = decode_record_components(packet, rec, groups_per_row, components, row_multiplier)\r
        for ci in range(components):\r
            start = ci * row_multiplier\r
            end = start + row_multiplier\r
            planes[ci].extend(row[:width] for row in comp_rows[start:end])\r
        final_states.append(states)\r
\r
    arrays = [np.asarray(rows, dtype=np.int32) for rows in planes]\r
    meta = {\r
        "packet": {"group": group, "index": index, "type": packet_type},\r
        "stream_index": stream_index,\r
        "shape": [int(arrays[0].shape[0]), int(arrays[0].shape[1])],\r
        "row_multiplier": row_multiplier,\r
        "components": components,\r
        "final_state": final_states[-1] if final_states else [],\r
    }\r
    return arrays, meta\r
\r
\r
def decode_type1_packet(arw: Path, group: int, index: int, out: Path | None = None, stream_index: int = 0) -> dict[str, object]:\r
    """Decode an independently parsed type-1 packet into int32 coefficient rows."""\r
\r
    packet, info = load_packet(arw, group, index, stream_index=stream_index)\r
    if info["type2"] != 1:\r
        raise ValueError(f"packet g{group}i{index} is type {info['type2']}, not type 1")\r
    width = infer_packet_width(group, info["type2"], int(info.get("work_width", info.get("raw_width", 7040))))\r
    groups_per_row = (width + 3) // 4\r
    rows: list[list[int]] = []\r
    final_states: list[dict[str, int]] = []\r
    for rec in info["records"]:\r
        if rec["byte_length"] <= 0:\r
            coeffs = [0] * (groups_per_row * 4)\r
            br = NativeBitReader.for_record(packet, rec["payload_offset"], rec["byte_length"])\r
        else:\r
            coeffs, _widths, br, _initial = decode_type1_row(packet, rec, groups_per_row)\r
        rows.append(coeffs[:width])\r
        final_states.append({"ptr": br.ptr, "bit": br.bit, "words_left": br.words_left, "status": br.status})\r
    if out:\r
        import numpy as np\r
\r
        arr = np.asarray(rows, dtype=np.int32)\r
        out.parent.mkdir(parents=True, exist_ok=True)\r
        arr.tofile(out)\r
    return {\r
        "packet": {"group": group, "index": index, "type": info["type2"]},\r
        "shape": [len(rows), width],\r
        "first_row_first16": rows[0][:16],\r
        "last_nonempty_row_first16": next((rows[i][:16] for i in range(len(rows) - 1, -1, -1) if any(rows[i])), []),\r
        "final_states_tail": final_states[-5:],\r
        "out": str(out) if out else "",\r
    }\r
\r
\r
def integrate_type1_coefficients(coeffs: "object", dc_offset: int) -> "object":\r
    """Apply the row postprocess after the entropy call at 0x1a7fb3.\r
\r
    For each row:\r
\r
      acc0 = int16(coeff[0]) * 2\r
      out[0] = acc0 >> 1\r
      acc_i = acc_{i-1} + int16(coeff[i]) * 2\r
      out[i] = acc_i >> 1\r
\r
    Applying this to group0/index0 and adding 2048 reproduces the native v4 c0\r
    lowpass plane.\r
    """\r
\r
    import numpy as np\r
\r
    c = np.asarray(coeffs, dtype=np.int32)\r
    out = np.empty_like(c, dtype=np.int32)\r
    signed = c.astype(np.int16).astype(np.int32)\r
    acc = signed[:, 0] * 2\r
    out[:, 0] = acc >> 1\r
    for x in range(1, c.shape[1]):\r
        acc = acc + signed[:, x] * 2\r
        out[:, x] = acc >> 1\r
    return (out + dc_offset).astype(np.int32)\r
\r
\r
def main() -> None:\r
    ap = argparse.ArgumentParser()\r
    ap.add_argument("arw", nargs="?", default="DSC00089.ARW")\r
    ap.add_argument("--group", type=int, default=0)\r
    ap.add_argument("--index", type=int, default=2)\r
    ap.add_argument("--row", type=int, default=280)\r
    ap.add_argument("--groups", type=int, default=16)\r
    ap.add_argument("--packet", action="store_true", help="decode the whole type-1 packet instead of one row")\r
    ap.add_argument("--components", action="store_true", help="decode a type-1/type-3 packet into component files")\r
    ap.add_argument("--out", default="")\r
    ns = ap.parse_args()\r
    if ns.components:\r
        result = decode_packet_components(Path(ns.arw), ns.group, ns.index, Path(ns.out) if ns.out else None)\r
    elif ns.packet:\r
        result = decode_type1_packet(Path(ns.arw), ns.group, ns.index, Path(ns.out) if ns.out else None)\r
    else:\r
        result = replay_row(Path(ns.arw), ns.group, ns.index, ns.row, ns.groups)\r
    text = json.dumps(result, indent=2)\r
    if ns.out and not ns.packet and not ns.components:\r
        Path(ns.out).write_text(text, encoding="utf-8")\r
    print(text)\r
\r
\r
if __name__ == "__main__":\r
    main()\r
`,Ar=`#!/usr/bin/env python3\r
"""Integer helpers from the LLVC3 reversing notes.\r
\r
Closer to a lab notebook than a polished codec module. The small functions\r
below are named around the traces they came from, then reused by the pure\r
decoder once a stage lines up.\r
"""\r
\r
from __future__ import annotations\r
\r
import argparse\r
import json\r
from pathlib import Path\r
\r
import numpy as np\r
\r
\r
INTERNAL_BIAS = 2048\r
\r
\r
def clamp_sample(x: np.ndarray, bits: int) -> np.ndarray:\r
    return np.clip(x, 0, (1 << bits) - 1)\r
\r
\r
def merge_average_detail(lo0: np.ndarray, lo1: np.ndarray, detail: np.ndarray, bits: int = 16) -> np.ndarray:\r
    """Inverse merge from Edit.exe RVA 0x1a9a40.\r
\r
    The scalar path is simple enough to leave here as a breadcrumb:\r
\r
        avg = (lo0 + lo1) >> 1\r
        sample = avg + 2 * signed_detail\r
        sample = clamp(sample, 0, (1 << bits) - 1)\r
\r
    Native output buffers are 16-bit, even when the TIFF tags describe 14-bit\r
    sensor samples, so the helper returns uint16 too.\r
    """\r
\r
    x = ((lo0.astype(np.int32) + lo1.astype(np.int32)) >> 1) + 2 * detail.astype(np.int32)\r
    return clamp_sample(x, bits).astype(np.uint16)\r
\r
\r
def add_detail(low: np.ndarray, detail: np.ndarray, bits: int = 16) -> np.ndarray:\r
    """Lowpass plus signed detail, from the 0x1aafd0 helper family."""\r
\r
    x = low.astype(np.int32) + detail.astype(np.int32)\r
    return clamp_sample(x, bits).astype(np.uint16)\r
\r
\r
def add_double_detail(low: np.ndarray, detail: np.ndarray, bits: int = 16) -> np.ndarray:\r
    """Lowpass plus 2*detail, seen around the 0x1ab2b0 helpers."""\r
\r
    x = low.astype(np.int32) + 2 * detail.astype(np.int32)\r
    return clamp_sample(x, bits).astype(np.uint16)\r
\r
\r
def sony_inv53_1d(low: np.ndarray, high: np.ndarray, axis: int) -> np.ndarray:\r
    """One-axis inverse 5/3 lifting in Sony's signed working domain."""\r
\r
    lo = np.asarray(low, dtype=np.int32)\r
    hi = np.asarray(high, dtype=np.int32)\r
    if axis == 0:\r
        hi_prev = np.vstack([hi[:1], hi[:-1]])\r
        lo2 = lo - ((hi_prev + hi + 2) >> 2)\r
        lo_next = np.vstack([lo2[1:], lo2[-1:]])\r
        hi2 = hi + ((lo2 + lo_next) >> 1)\r
        out = np.empty((lo.shape[0] * 2, lo.shape[1]), dtype=np.int32)\r
        out[0::2] = lo2\r
        out[1::2] = hi2\r
        return out\r
    if axis == 1:\r
        hi_prev = np.concatenate([hi[:, :1], hi[:, :-1]], axis=1)\r
        lo2 = lo - ((hi_prev + hi + 2) >> 2)\r
        lo_next = np.concatenate([lo2[:, 1:], lo2[:, -1:]], axis=1)\r
        hi2 = hi + ((lo2 + lo_next) >> 1)\r
        out = np.empty((lo.shape[0], lo.shape[1] * 2), dtype=np.int32)\r
        out[:, 0::2] = lo2\r
        out[:, 1::2] = hi2\r
        return out\r
    raise ValueError("axis must be 0 or 1")\r
\r
\r
def sony_inv53_1d_high_leading(low: np.ndarray, high: np.ndarray) -> np.ndarray:\r
    """Vertical inverse 5/3 where Sony's guard line makes the high row lead."""\r
\r
    lo = np.asarray(low, dtype=np.int32)\r
    hi = np.asarray(high, dtype=np.int32)\r
    if hi.shape[1] != lo.shape[1]:\r
        raise ValueError(f"unexpected high-leading shapes: low={lo.shape}, high={hi.shape}")\r
\r
    if hi.shape[0] == lo.shape[0] + 1:\r
        lo2 = lo - ((hi[:-1] + hi[1:] + 2) >> 2)\r
        hi2 = np.empty_like(hi)\r
        hi2[0] = hi[0] + lo2[0]\r
        if lo2.shape[0] > 1:\r
            hi2[1:-1] = hi[1:-1] + ((lo2[:-1] + lo2[1:]) >> 1)\r
        hi2[-1] = hi[-1] + lo2[-1]\r
        out = np.empty((lo.shape[0] * 2 + 1, lo.shape[1]), dtype=np.int32)\r
    elif hi.shape[0] == lo.shape[0]:\r
        hi_next = np.vstack([hi[1:], hi[-1:]])\r
        lo2 = lo - ((hi + hi_next + 2) >> 2)\r
        hi2 = np.empty_like(hi)\r
        hi2[0] = hi[0] + lo2[0]\r
        if lo2.shape[0] > 1:\r
            hi2[1:] = hi[1:] + ((lo2[:-1] + lo2[1:]) >> 1)\r
        out = np.empty((lo.shape[0] * 2, lo.shape[1]), dtype=np.int32)\r
    else:\r
        raise ValueError(f"unexpected high-leading row counts: low={lo.shape}, high={hi.shape}")\r
\r
    out[0::2] = hi2\r
    out[1::2] = lo2\r
    return out\r
\r
\r
def llvc3_edge_detail(x: np.ndarray, edge_mode: str = "even") -> np.ndarray:\r
    """Expand Sony's edge-only HH detail row."""\r
\r
    xi = x.astype(np.int32)\r
    signed_half_step = np.where(xi > 0, 1, np.where(xi < 0, -1, 0))\r
    if edge_mode == "even":\r
        mask = ((xi & 1) == 0) & (xi != 0)\r
    elif edge_mode == "odd":\r
        mask = (xi & 1) != 0\r
    else:\r
        raise ValueError(f"unknown LLVC3 edge_mode {edge_mode!r}")\r
    return (2 * xi + np.where(mask, signed_half_step, 0)).astype(np.int32)\r
\r
\r
def synthesize_llvc3_level(ll: np.ndarray, sub0: np.ndarray, sub1: np.ndarray, sub2: np.ndarray) -> np.ndarray:\r
    """Synthesize one LLVC3 scale from LL plus three detail subbands.\r
\r
    Packet component mapping, as verified against Imaging Edge for this ARW6\r
    sample:\r
\r
    * sub0: horizontal detail for the low vertical branch (HL)\r
    * sub1: vertical detail for the low horizontal branch (LH)\r
    * sub2: diagonal detail (HH)\r
\r
    The annoying bit is the extra flush row. At the bottom edge Sony feeds it\r
    across the branches: LH's final row comes from sub0[-1], HH's final row from\r
    sub1[-1].\r
    """\r
\r
    ll_i = np.asarray(ll, dtype=np.int32)\r
    h = ll_i.shape[0]\r
    if ll_i.shape[1] != sub0.shape[1] or sub0.shape != sub1.shape or sub1.shape != sub2.shape:\r
        raise ValueError(f"unexpected subband shapes: ll={ll_i.shape}, sub0={sub0.shape}, sub1={sub1.shape}, sub2={sub2.shape}")\r
    if sub0.shape[0] < h + 1:\r
        raise ValueError(f"subbands need one flush row: ll={ll_i.shape}, sub={sub0.shape}")\r
\r
    lh = np.empty((h, ll_i.shape[1]), dtype=np.int32)\r
    lh[:-1] = sub1[1:h]\r
    lh[-1] = sub0[h]\r
\r
    hh = np.empty((h, ll_i.shape[1]), dtype=np.int32)\r
    hh[:-1] = sub2[1:h]\r
    hh[-1] = sub1[h]\r
\r
    low_horizontal = sony_inv53_1d(ll_i, lh, axis=0)\r
    high_horizontal = sony_inv53_1d(sub0[:h], hh, axis=0)\r
    return sony_inv53_1d(low_horizontal, high_horizontal, axis=1)\r
\r
\r
def trunc_div2(x: np.ndarray) -> np.ndarray:\r
    """Integer division by two with C/C++ truncation toward zero."""\r
\r
    a = np.asarray(x, dtype=np.int32)\r
    return np.where(a >= 0, a // 2, -((-a) // 2)).astype(np.int32)\r
\r
\r
def synthesize_llvc3_level_stride(\r
    ll: np.ndarray,\r
    sub0: np.ndarray,\r
    sub1: np.ndarray,\r
    sub2: np.ndarray,\r
    edge_rows: int,\r
    bottom_hh_extra: np.ndarray | None = None,\r
    edge_mode: str = "even",\r
) -> np.ndarray:\r
    """Same synthesis, with the larger line-flush padding used above group 1.\r
\r
    \`edge_rows\` is the row multiplier minus one for the higher scales:\r
    group1 -> 0, group2 -> 1, group3 -> 2.  The group1 case falls back to the\r
    smaller helper above.\r
    """\r
\r
    if edge_rows == 0:\r
        return synthesize_llvc3_level(ll, sub0, sub1, sub2)\r
\r
    ll_i = np.asarray(ll, dtype=np.int32)\r
    h, w = ll_i.shape\r
    if sub0.shape != sub1.shape or sub1.shape != sub2.shape or sub0.shape[1] != w:\r
        raise ValueError(f"unexpected subband shapes: ll={ll_i.shape}, sub0={sub0.shape}, sub1={sub1.shape}, sub2={sub2.shape}")\r
    if sub0.shape[0] < h + edge_rows * 2:\r
        raise ValueError(f"not enough line-flush rows: ll={ll_i.shape}, sub={sub0.shape}, edge_rows={edge_rows}")\r
\r
    hl = np.empty((h, w), dtype=np.int32)\r
    hl[:edge_rows] = sub0[:edge_rows]\r
    hl[edge_rows : h - edge_rows] = sub0[2 * edge_rows : h]\r
    hl[h - edge_rows :] = sub0[h : h + edge_rows]\r
\r
    lh = np.empty((h, w), dtype=np.int32)\r
    lh[:edge_rows] = sub0[edge_rows : 2 * edge_rows]\r
    lh[edge_rows : h - edge_rows] = sub1[2 * edge_rows : h]\r
    lh[h - edge_rows :] = sub0[h + edge_rows : h + 2 * edge_rows]\r
\r
    hh = np.empty((h, w), dtype=np.int32)\r
\r
    hh[:edge_rows] = llvc3_edge_detail(sub1[:edge_rows], edge_mode)\r
    hh[edge_rows : h - edge_rows] = sub2[np.arange(edge_rows, h - edge_rows) + edge_rows]\r
    bottom_hh = llvc3_edge_detail(sub1[h : h + edge_rows], edge_mode)\r
    if bottom_hh_extra is not None:\r
        extra = np.asarray(bottom_hh_extra, dtype=np.int32)\r
        if extra.shape != bottom_hh.shape:\r
            raise ValueError(f"bottom_hh_extra shape {extra.shape} != bottom edge {bottom_hh.shape}")\r
        bottom_hh = bottom_hh + extra\r
    hh[h - edge_rows :] = bottom_hh\r
\r
    low_horizontal = sony_inv53_1d(ll_i, lh, axis=0)\r
    high_horizontal = sony_inv53_1d(hl, hh, axis=0)\r
    return sony_inv53_1d(low_horizontal, high_horizontal, axis=1)\r
\r
\r
def synthesize_llvc3_guard_group1(ll: np.ndarray, sub0: np.ndarray, sub1: np.ndarray, sub2: np.ndarray) -> np.ndarray:\r
    """Guard-row group 1 synthesis used by non-16-aligned LLVC heights."""\r
\r
    ll_i = np.asarray(ll, dtype=np.int32)\r
    h, w = ll_i.shape\r
    if sub0.shape != sub1.shape or sub1.shape != sub2.shape or sub0.shape[1] != w:\r
        raise ValueError(f"unexpected subband shapes: ll={ll_i.shape}, sub0={sub0.shape}, sub1={sub1.shape}, sub2={sub2.shape}")\r
    if sub0.shape[0] < h + 2:\r
        raise ValueError(f"not enough guarded group1 rows: ll={ll_i.shape}, sub={sub0.shape}")\r
\r
    lh = np.empty((h + 1, w), dtype=np.int32)\r
    lh[:-1] = sub1[1 : 1 + h]\r
    lh[-1] = sub0[h + 1]\r
\r
    hh = np.empty((h + 1, w), dtype=np.int32)\r
    hh[:-1] = sub2[1 : 1 + h]\r
    hh[-1] = sub1[h + 1]\r
\r
    low_horizontal = sony_inv53_1d_high_leading(ll_i, lh)\r
    high_horizontal = sony_inv53_1d_high_leading(sub0[1 : 1 + h], hh)\r
    return sony_inv53_1d(low_horizontal, high_horizontal, axis=1)\r
\r
\r
def synthesize_llvc3_guard_group2(\r
    ll: np.ndarray, sub0: np.ndarray, sub1: np.ndarray, sub2: np.ndarray, edge_mode: str = "even"\r
) -> np.ndarray:\r
    """Guard-row group 2 synthesis for cropped-height ARW6 tiles."""\r
\r
    ll_i = np.asarray(ll, dtype=np.int32)\r
    h, w = ll_i.shape\r
    if sub0.shape != sub1.shape or sub1.shape != sub2.shape or sub0.shape[1] != w:\r
        raise ValueError(f"unexpected subband shapes: ll={ll_i.shape}, sub0={sub0.shape}, sub1={sub1.shape}, sub2={sub2.shape}")\r
    if sub0.shape[0] < h + 3:\r
        raise ValueError(f"not enough guarded group2 rows: ll={ll_i.shape}, sub={sub0.shape}")\r
\r
    hl = sub0[2 : 2 + h]\r
    lh = np.empty((h, w), dtype=np.int32)\r
    lh[0] = sub0[0]\r
    lh[1:] = sub1[2 : 1 + h]\r
\r
    hh = np.empty((h, w), dtype=np.int32)\r
    hh[0] = llvc3_edge_detail(sub0[1:2], edge_mode)[0]\r
    hh[1:] = sub2[2 : 1 + h]\r
\r
    low_horizontal = sony_inv53_1d_high_leading(ll_i, lh)\r
    high_horizontal = sony_inv53_1d_high_leading(hl, hh)\r
    return sony_inv53_1d(low_horizontal, high_horizontal, axis=1)\r
\r
\r
def synthesize_llvc3_guard_group3(\r
    ll: np.ndarray, sub0: np.ndarray, sub1: np.ndarray, sub2: np.ndarray, edge_mode: str = "even"\r
) -> np.ndarray:\r
    """Guard-row group 3 synthesis for cropped-height ARW6 tiles."""\r
\r
    ll_i = np.asarray(ll, dtype=np.int32)\r
    h, w = ll_i.shape\r
    if sub0.shape != sub1.shape or sub1.shape != sub2.shape or sub0.shape[1] != w:\r
        raise ValueError(f"unexpected subband shapes: ll={ll_i.shape}, sub0={sub0.shape}, sub1={sub1.shape}, sub2={sub2.shape}")\r
    if sub0.shape[0] < h + 5:\r
        raise ValueError(f"not enough guarded group3 rows: ll={ll_i.shape}, sub={sub0.shape}")\r
\r
    hl = np.empty((h, w), dtype=np.int32)\r
    hl[0] = sub0[0]\r
    hl[1:] = sub0[4 : 3 + h]\r
\r
    lh = np.empty((h, w), dtype=np.int32)\r
    lh[0] = sub0[1]\r
    lh[1:-1] = sub1[4 : 2 + h]\r
    lh[-1] = sub0[h + 3]\r
\r
    hh = np.empty((h, w), dtype=np.int32)\r
    hh[0] = llvc3_edge_detail(sub0[2:3], edge_mode)[0]\r
    hh[1:-1] = sub2[4 : 2 + h]\r
    hh[-1] = llvc3_edge_detail(sub0[h + 4 : h + 5], edge_mode)[0]\r
\r
    low_horizontal = sony_inv53_1d(ll_i, lh, axis=0)\r
    high_horizontal = sony_inv53_1d(hl, hh, axis=0)\r
    return sony_inv53_1d(low_horizontal, high_horizontal, axis=1)\r
\r
\r
def synthesize_llvc3_final_green(ll: np.ndarray, detail: np.ndarray, top_rows: int = 4) -> np.ndarray:\r
    """Final CFA-green reconstruction, from the 0x1ab570 path.\r
\r
    Not the same 2-D 5/3 inverse used by groups 1..3. It expands the half-width\r
    green lowpass into both RGGB green sites. The row offsets look odd because\r
    the native line buffer keeps four guard rows at the top and a few latency\r
    rows at the bottom.\r
    """\r
\r
    ll_i = np.asarray(ll, dtype=np.int32)\r
    det = np.asarray(detail, dtype=np.int32)\r
    h, w = ll_i.shape\r
    if not 0 <= top_rows <= 8:\r
        raise ValueError(f"unexpected final green top row count {top_rows}")\r
    if det.shape[1] != w or det.shape[0] < 8 + max(0, h - top_rows):\r
        raise ValueError(f"unexpected final green shapes: ll={ll_i.shape}, detail={det.shape}")\r
\r
    selected = np.empty((h, w), dtype=np.int32)\r
    top = min(top_rows, h)\r
    selected[:top] = det[:top]\r
    if h > top:\r
        selected[top:] = det[8 : 8 + (h - top)]\r
\r
    odd_green = np.empty((h, w), dtype=np.int32)\r
    for y in range(h):\r
        cur = selected[y]\r
        prev = selected[y - 1] if y > 0 else cur\r
        pred = np.empty(w, dtype=np.int32)\r
        pred[:-1] = (cur[1:] + prev[:-1] + cur[:-1] + prev[1:]) >> 2\r
        pred[-1] = ((prev[-1] + cur[-1]) * 2) >> 2\r
        odd_green[y] = ((2 * ll_i[y] - pred) >> 1).astype(np.int32)\r
\r
    even_green = np.empty((h, w), dtype=np.int32)\r
    for y in range(h):\r
        cur = odd_green[y]\r
        nxt = odd_green[y + 1] if y + 1 < h else cur\r
        even_green[y, 0] = selected[y, 0] + (((cur[0] + nxt[0]) * 2) >> 2)\r
        if w > 1:\r
            even_green[y, 1:] = selected[y, 1:] + ((cur[:-1] + nxt[:-1] + nxt[1:] + cur[1:]) >> 2)\r
\r
    out = np.empty((h, w * 2), dtype=np.int32)\r
    out[:, 0::2] = even_green\r
    out[:, 1::2] = odd_green\r
    return out\r
\r
\r
def finalize_llvc3_color_planes(\r
    v1_green: np.ndarray, v1_red: np.ndarray, v1_blue: np.ndarray, full_green: np.ndarray\r
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:\r
    """Apply the final color-plane relation for the decoded output.\r
\r
    Red and blue keep the group-3 residuals. Only the green predictor is swapped\r
    for the final CFA-green average:\r
\r
        residual = (v1_color - v1_green) / 2\r
        v0_color = average(clamp12(final_green_pair)) + 2 * residual\r
\r
    Sony clamps the two final green sites to the 12-bit code domain before using\r
    them as the red/blue predictor. Without that, highlight green overshoot leaks\r
    into R/B as one to three LUT code steps.\r
    """\r
\r
    g = np.asarray(v1_green, dtype=np.int32)\r
    r = np.asarray(v1_red, dtype=np.int32)\r
    b = np.asarray(v1_blue, dtype=np.int32)\r
    fg = np.asarray(full_green, dtype=np.int32)\r
    if fg.shape != (g.shape[0], g.shape[1] * 2):\r
        raise ValueError(f"unexpected final green shape: v1={g.shape}, full={fg.shape}")\r
    if r.shape != g.shape or b.shape != g.shape:\r
        raise ValueError(f"unexpected v1 color shapes: green={g.shape}, red={r.shape}, blue={b.shape}")\r
\r
    fg_pred = np.clip(fg + INTERNAL_BIAS, 0, 4095).astype(np.int32) - INTERNAL_BIAS\r
    gavg = (fg_pred[:, 0::2] + fg_pred[:, 1::2]) >> 1\r
    red_residual = (r - g) // 2\r
    blue_residual = (b - g) // 2\r
    return fg, gavg + 2 * red_residual, gavg + 2 * blue_residual\r
\r
\r
def signed_to_sample(x: np.ndarray, bits: int = 16, bias: int = INTERNAL_BIAS) -> np.ndarray:\r
    """Convert LLVC3 signed internal rows to Sony's unsigned output samples."""\r
\r
    return clamp_sample(np.asarray(x, dtype=np.int32) + bias, bits).astype(np.uint16)\r
\r
\r
def apply_sample_lut(code_samples: np.ndarray, lut: np.ndarray) -> np.ndarray:\r
    """Map unsigned LLVC3 code-domain samples through a Sony sample LUT."""\r
\r
    table = np.asarray(lut, dtype=np.uint16).reshape(-1)\r
    if table.size == 0:\r
        raise ValueError("sample LUT is empty")\r
    code = np.clip(np.asarray(code_samples, dtype=np.int32), 0, table.size - 1)\r
    return table[code].astype(np.uint16)\r
\r
\r
def clamp_signed_to_code_range(x: np.ndarray, max_code: int = 4095, bias: int = INTERNAL_BIAS) -> np.ndarray:\r
    """Clamp signed LLVC3 rows to Sony's 12-bit code range, then return signed rows."""\r
\r
    return np.clip(np.asarray(x, dtype=np.int32) + bias, 0, max_code).astype(np.int32) - bias\r
\r
\r
def lifting_predict_detail(detail: np.ndarray, a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray) -> np.ndarray:\r
    """Prediction/update kernel from the 0x1ab570 scalar path.\r
\r
    The four neighbor names are still placeholder-ish, but the integer\r
    operation itself is clear in the trace:\r
\r
        pred = (a + b + c + d) >> 2\r
        out = (2 * detail - pred) >> 1\r
    """\r
\r
    pred = (a.astype(np.int32) + b.astype(np.int32) + c.astype(np.int32) + d.astype(np.int32)) >> 2\r
    return ((2 * detail.astype(np.int32) - pred) >> 1).astype(np.int32)\r
\r
\r
def recombine_rggb(c0: np.ndarray, c1: np.ndarray, c2: np.ndarray) -> np.ndarray:\r
    """Recombine decoded LLVC planes into the TIFF-declared RGGB Bayer mosaic."""\r
\r
    half_h, width = c0.shape\r
    if c1.shape != (half_h, width // 2) or c2.shape != (half_h, width // 2):\r
        raise ValueError(f"unexpected plane shapes: {c0.shape}, {c1.shape}, {c2.shape}")\r
    out = np.empty((half_h * 2, width), dtype=np.uint16)\r
    out[0::2, 0::2] = c1\r
    out[0::2, 1::2] = c0[:, 1::2]\r
    out[1::2, 0::2] = c0[:, 0::2]\r
    out[1::2, 1::2] = c2\r
    return out\r
\r
\r
def selftest() -> dict[str, object]:\r
    lo0 = np.array([1000, 1002, 10, 65530], dtype=np.uint16)\r
    lo1 = np.array([1002, 1004, 10, 65530], dtype=np.uint16)\r
    detail = np.array([0, 1, -20, 20], dtype=np.int32)\r
    merged = merge_average_detail(lo0, lo1, detail, bits=16)\r
    added = add_detail(lo0, detail, bits=16)\r
    doubled = add_double_detail(lo0, detail, bits=16)\r
    pred = lifting_predict_detail(\r
        detail,\r
        np.array([4, 8, 12, 16]),\r
        np.array([4, 8, 12, 16]),\r
        np.array([4, 8, 12, 16]),\r
        np.array([4, 8, 12, 16]),\r
    )\r
    return {\r
        "merge_average_detail": merged.tolist(),\r
        "add_detail": added.tolist(),\r
        "add_double_detail": doubled.tolist(),\r
        "lifting_predict_detail": pred.tolist(),\r
    }\r
\r
\r
def main() -> None:\r
    ap = argparse.ArgumentParser()\r
    ap.add_argument("--selftest", action="store_true")\r
    ap.add_argument("--out", default="")\r
    ns = ap.parse_args()\r
    result = selftest()\r
    text = json.dumps(result, indent=2)\r
    if ns.out:\r
        Path(ns.out).write_text(text, encoding="utf-8")\r
    print(text)\r
\r
\r
if __name__ == "__main__":\r
    main()\r
`;const Lt=512,Ir=kr["Sony ILCE-7RM5"].split(/\s+/).map(Number).filter(Number.isFinite);let lt=null;function Tr(n){if(n.byteLength<8)return null;const e=n.getUint16(0,!1);return e===18761?!0:e===19789?!1:null}function rr(n,e,t){const r=Math.min(n.length,e+t);let i="";for(let s=e;s<r;s++){const a=n[s];if(a===0)break;i+=String.fromCharCode(a)}return i.trim()}function K(n,e,t,r,i){const s=t===1||t===2||t===7?1:t===3||t===8?2:t===4||t===9?4:0;if(!s)return[];const a=s*r,o=a<=4?i:n.getUint32(i,e);if(o<0||o+a>n.byteLength)return[];const l=[];for(let h=0;h<r;h++){const c=o+h*s;t===1||t===2||t===7?l.push(n.getUint8(c)):t===3?l.push(n.getUint16(c,e)):t===8?l.push(n.getInt16(c,e)):t===4?l.push(n.getUint32(c,e)):t===9&&l.push(n.getInt32(c,e))}return l}function Ot(n,e,t,r,i,s){if(r!==2||i<=0)return"";const a=i<=4?s:e.getUint32(s,t);return a<0||a>=n.length?"":rr(n,a,i)}function nr(n){const e=new Uint8Array(n),t=new DataView(n),r=Tr(t);if(r===null||t.getUint16(2,r)!==42)return null;const s=d=>t.getUint16(d,r),a=d=>t.getUint32(d,r),o=[a(4)],l=new Set;let h="",c="";for(;o.length;){const d=o.pop();if(l.has(d)||d<=0||d+2>t.byteLength)continue;l.add(d);const u=s(d);if(d+2+u*12+4>t.byteLength)continue;const p=new Map;for(let T=0;T<u;T++){const P=d+2+T*12,O=s(P),E=s(P+2),U=a(P+4),k=P+8;p.set(O,{type:E,count:U,valueOffset:k})}const _=p.get(271),g=p.get(272);_&&!h&&(h=Ot(e,t,r,_.type,_.count,_.valueOffset)),g&&!c&&(c=Ot(e,t,r,g.type,g.count,g.valueOffset));const y=p.get(330);if(y){const T=K(t,r,y.type,y.count,y.valueOffset);for(const P of T)o.push(P)}const b=p.get(259),v=p.get(262);if(b&&v){const T=K(t,r,b.type,b.count,b.valueOffset)[0],P=K(t,r,v.type,v.count,v.valueOffset)[0];if(T===32766&&P===32803){const O=K(t,r,p.get(256).type,p.get(256).count,p.get(256).valueOffset)[0],E=K(t,r,p.get(257).type,p.get(257).count,p.get(257).valueOffset)[0],U=K(t,r,p.get(258).type,p.get(258).count,p.get(258).valueOffset)[0],k=K(t,r,p.get(273).type,p.get(273).count,p.get(273).valueOffset)[0],W=K(t,r,p.get(279).type,p.get(279).count,p.get(279).valueOffset)[0],q=p.get(33422)?K(t,r,p.get(33422).type,p.get(33422).count,p.get(33422).valueOffset):[0,1,1,2];p.get(29456)&&K(t,r,p.get(29456).type,p.get(29456).count,p.get(29456).valueOffset);const X=p.get(50717)?K(t,r,p.get(50717).type,p.get(50717).count,p.get(50717).valueOffset)[0]:16383,x=p.get(50719)?K(t,r,p.get(50719).type,p.get(50719).count,p.get(50719).valueOffset):[],B=p.get(50720)?K(t,r,p.get(50720).type,p.get(50720).count,p.get(50720).valueOffset):[];if(k+Lt+16>e.length||k+W>e.length)return null;const f=k+Lt,C=rr(e,f,4),D=e[f+8]<<8|e[f+9],S=e[f+10]<<8|e[f+11],w=e[f+12]<<8|e[f+13],R=e[f+14]<<8|e[f+15],J=w>>4&63,ue=R>>13,re=R>>10&3,Q=S*2,V=W>=4?(e[k]|e[k+1]<<8|e[k+2]<<16|e[k+3]<<24)>>>0:0,rt=D===O&&Q===E;let ve=!1;if(V>=1&&V<=16&&W>=8+V*24){const ne=new Map,se=new Map;let Te=!0;for(let Re=0;Re<V;Re++){const ae=k+8+Re*24,fe=t.getUint32(ae+8,!0),ge=t.getUint32(ae+12,!0),st=t.getUint32(ae+16,!0),Ee=t.getUint32(ae+20,!0);if(!st||!Ee||fe+st>O||ge+Ee>E){Te=!1;break}const Tt=ne.get(ge);if(Tt!==void 0&&Tt!==Ee){Te=!1;break}ne.set(ge,Ee),se.set(ge,(se.get(ge)||0)+st)}if(Te){const Re=Array.from(ne.keys()).sort((fe,ge)=>fe-ge);let ae=0;for(const fe of Re){if(fe!==ae||se.get(fe)!==O){Te=!1;break}ae+=ne.get(fe)}ve=Te&&ae===E}}const Me=V>=1&&V<=16&&D>0&&Q>0&&O%D===0&&E%Q===0&&V===O/D*(E/Q);if(C!=="A000"&&C!=="0000"||!rt&&!Me&&!ve||J!==16||ue!==3||re!==3)return null;const nt=[1024,1024,1024,1024],it=q.length>=4?q.slice(0,4).map(ne=>ne===0?"R":ne===2?"B":ne===1?"G":"?").join(""):"";return{width:O,height:E,bitsPerSample:U,compression:T,photometric:P,blackLevel:nt,whiteLevel:Number(X||16383),cfaPattern:it,defaultCropOrigin:x.length>=2?[Number(x[0]),Number(x[1])]:void 0,defaultCropSize:B.length>=2?[Number(B[0]),Number(B[1])]:void 0,make:h||"SONY",model:c||"ILCE-7M5"}}}const A=a(d+2+u*12);A&&o.push(A)}return null}async function Rr(n){return lt||(lt=(async()=>{const e=await Cr();return e.__jtrSonyCrawHqDecoderReady||(await e.FS.mkdirTree("/sony_craw_hq"),await e.FS.writeFile("/sony_craw_hq/llvc3_bitstream_probe.py",vr),await e.FS.writeFile("/sony_craw_hq/llvc3_entropy.py",Pr),await e.FS.writeFile("/sony_craw_hq/llvc3_math.py",Ar),await e.runPythonAsync(`
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
`),e.__jtrSonyCrawHqDecoderReady=!0),e})()),lt}function Lr(n){return nr(n)}function Or(n){return n==="ILCE-7M5"?Ir:null}function ir(n,e,t,r){var l;if(n.length!==e.width*e.height)throw new Error(`Sony cRAW HQ decoded size mismatch: got ${n.length}, expected ${e.width*e.height}`);const i=e.model||"ILCE-7M5",s=i.startsWith("Sony ")?i:`Sony ${i}`,a=Or(i),o={...r||{},make:e.make||(r==null?void 0:r.camera_make)||"SONY",model:i,camera_make:e.make||(r==null?void 0:r.camera_make)||"SONY",camera_model:i,UniqueCameraModel:s,sourceFormat:t==="joraw2-wasm"?"Sony cRAW HQ / LLVC3 (JoRaw2 WASM)":"Sony cRAW HQ / LLVC3",sonyCrawHq:{...e,decodeBackend:t},color_desc:e.cfaPattern,black_level_per_channel:e.blackLevel,white_level:e.whiteLevel,color_matrix:a&&a.length===9?a:void 0,idata:{filters:2492765332,colors:3},color_data:{...(r==null?void 0:r.color_data)||{},black:1024,cblack_rawpy_style:e.blackLevel,dng_levels:{...((l=r==null?void 0:r.color_data)==null?void 0:l.dng_levels)||{},dng_cblack:e.blackLevel,dng_whitelevel:e.whiteLevel}}};return{data:n,width:e.width,height:e.height,bayerPattern:e.cfaPattern,blackLevels:e.blackLevel,whiteLevel:e.whiteLevel,metadata:o,isThreePlane:!1,isXTrans:!1}}async function Dr(n,e,t){const r=typeof performance<"u"?performance.now():Date.now(),i=await bt(),s=typeof performance<"u"?performance.now():Date.now(),a=new i;try{const o=new Uint8Array(n);await a.open(o,{});const l=typeof performance<"u"?performance.now():Date.now();let h=null;try{h=await a.metadata(!0)}catch(y){console.warn("[Sony cRAW HQ] fast WASM metadata read failed",y)}const c=typeof performance<"u"?performance.now():Date.now(),d=a.getRawImage(),u=typeof performance<"u"?performance.now():Date.now();if(!d||!d.data)throw new Error("Sony cRAW HQ LibRaw WASM returned no raw image");const p=d.data instanceof Uint16Array?d.data:new Uint16Array(d.data.buffer,d.data.byteOffset||0,d.data.byteLength/2),_=typeof performance<"u"?performance.now():Date.now();if(d.width!==e.width||d.height!==e.height)throw new Error(`Sony cRAW HQ LibRaw WASM dimensions mismatch: got ${d.width}x${d.height}, expected ${e.width}x${e.height}`);const g=typeof performance<"u"?performance.now():Date.now();return console.info("[Sony cRAW HQ] fast decode timings",{width:e.width,height:e.height,backend:"joraw2-wasm",wasmReadyMs:Math.round(s-r),openMs:Math.round(l-s),metadataMs:Math.round(c-l),unpackMs:Math.round(u-c),copyMs:Math.round(_-u),totalMs:Math.round(g-r)}),{rawImageData:ir(p,e,"joraw2-wasm",h),info:e}}finally{typeof a.delete=="function"?a.delete():typeof a.close=="function"&&a.close()}}async function Fr(n,e,t){const r=typeof performance<"u"?performance.now():Date.now(),i=await Rr(),s=typeof performance<"u"?performance.now():Date.now(),a=new Uint8Array(n),o=await fetch(new URL("/assets/sony_llvc3_static_lut4096_padded_u16-FsVBk-IV.bin",import.meta.url));if(!o.ok)throw new Error(`Failed to load Sony LLVC3 sample LUT: HTTP ${o.status}`);const l=new Uint8Array(await o.arrayBuffer()),h=typeof performance<"u"?performance.now():Date.now();i.globals.set("jtr_sony_arw_bytes",a),i.globals.set("jtr_sony_lut_bytes",l);const c=await i.runPythonAsync("jtr_decode_sony_craw_hq(jtr_sony_arw_bytes.to_py(), jtr_sony_lut_bytes.to_py())"),d=typeof performance<"u"?performance.now():Date.now(),u=c.toJs();typeof c.destroy=="function"&&c.destroy(),i.globals.delete("jtr_sony_arw_bytes"),i.globals.delete("jtr_sony_lut_bytes");const p=new Uint8Array(u.byteLength);p.set(u);const _=new Uint16Array(p.buffer),g=typeof performance<"u"?performance.now():Date.now(),y=typeof performance<"u"?performance.now():Date.now();return console.info("[Sony cRAW HQ] decode timings",{width:e.width,height:e.height,backend:"pyodide",pyodideReadyMs:Math.round(s-r),lutLoadMs:Math.round(h-s),llvc3DecodeMs:Math.round(d-h),copyMs:Math.round(g-d),totalMs:Math.round(y-r)}),{rawImageData:ir(_,e,"pyodide"),info:e}}async function Mr(n,e){const t=nr(n);if(!t)return null;try{return await Dr(n,t,e)}catch(r){return console.warn("[Sony cRAW HQ] fast WASM decode failed; falling back to Pyodide",r),Fr(n,t)}}async function Er(n,e){return Mr(n,e)}const Ur=["RGGB","BGGR","GRBG","GBRG"],Nr=new Set(Ur);function Br(n){if(!n||typeof n!="object"||Array.isArray(n)||ArrayBuffer.isView(n))return n;const e=n;return e.value??e.values??e.description??n}function Ue(n){const e=Br(n);let t="";if(typeof e=="string")t=e.toUpperCase().replace(/[^RGB012]/g,"");else if(typeof e=="number")t=String(e);else if(Array.isArray(e)||ArrayBuffer.isView(e))t=Array.from(e).map(String).join("");else return null;return/^[012]{4}$/.test(t)&&(t=t.replace(/0/g,"R").replace(/1/g,"G").replace(/2/g,"B")),Nr.has(t)?t:null}function Dt(n){var e;return n?[n.cfa_pattern,n.cfaPattern,n.BayerPattern,n.CFAPattern2,n.CFAPattern,(e=n.idata)==null?void 0:e.cfa_pattern]:[]}function Vr(n,e,t){const r=Ue(n.bayerPattern),i=n.bayerPatternSource==="manual";if(r&&!i)return{pattern:r,source:n.bayerPatternSource||"decoder"};if(!i)for(const a of Dt(n.metadata)){const o=Ue(a);if(o)return{pattern:o,source:"libraw-metadata"}}for(const a of Dt(e)){const o=Ue(a);if(o)return{pattern:o,source:"metadata"}}const s=Ue(t)||(i?r:null);return s?{pattern:s,source:"manual"}:{pattern:null,source:null}}function zr(n,e){return n.bayerPattern=e.pattern||"",n.bayerPatternSource=e.source||void 0,e.pattern&&e.source!=="manual"&&(n.metadata={...n.metadata||{},cfa_pattern:e.pattern,cfaPattern:e.pattern}),e}const Z=n=>{const e=Number(n);return Number.isFinite(e)?Math.max(0,e):0};function me(n){if(!n||typeof n.length!="number")return null;const e=Array.from(n);return e.length<4?null:[Z(e[0]),Z(e[1]),Z(e[2]),Z(e[3])]}function Ft(n){if(!n||n.source!=="libraw")return null;const e=me(n.channelOffsets),t=me(n.channelLevels),r=me(n.siteColorIndices),i=me(n.siteBaseLevels),s=me(n.siteLevels);if(!e||!t||!r||!i||!s)return null;const a=Math.max(0,Math.floor(Z(n.repeatRows))),o=Math.max(0,Math.floor(Z(n.repeatCols))),l=a*o,h=n.repeatValues&&typeof n.repeatValues.length=="number"?Array.from(n.repeatValues,Z):[],c=l>0&&l<=4098&&h.length>=l;return{source:"libraw",common:Z(n.common),channelOffsets:e,channelLevels:t,siteColorIndices:r,siteBaseLevels:i,siteLevels:s,repeatRows:c?a:0,repeatCols:c?o:0,repeatOriginY:Math.max(0,Math.floor(Z(n.repeatOriginY))),repeatOriginX:Math.max(0,Math.floor(Z(n.repeatOriginX))),repeatValues:c?h.slice(0,l):[]}}var $e=typeof self<"u"?self:global;const Fe=typeof navigator<"u",Gr=Fe&&typeof HTMLImageElement>"u",je=!(typeof global>"u"||typeof process>"u"||!process.versions||!process.versions.node),Je=$e.Buffer,Ne=$e.BigInt,Qe=!!Je,Hr=n=>n;function We(n,e=Hr){if(je)try{return typeof require=="function"?Promise.resolve(e(require(n))):import(n).then(e)}catch{console.warn(`Couldn't load ${n}`)}}let xt=$e.fetch;const jr=n=>xt=n;if(!$e.fetch){const n=We("http",r=>r),e=We("https",r=>r),t=(r,{headers:i}={})=>new Promise(async(s,a)=>{let{port:o,hostname:l,pathname:h,protocol:c,search:d}=new URL(r);const u={method:"GET",hostname:l,path:encodeURI(h)+d,headers:i};o!==""&&(u.port=Number(o));const p=(c==="https:"?await e:await n).request(u,_=>{if(_.statusCode===301||_.statusCode===302){let g=new URL(_.headers.location,r).toString();return t(g,{headers:i}).then(s).catch(a)}s({status:_.statusCode,arrayBuffer:()=>new Promise(g=>{let y=[];_.on("data",b=>y.push(b)),_.on("end",()=>g(Buffer.concat(y)))})})});p.on("error",a),p.end()});jr(t)}function m(n,e,t){return e in n?Object.defineProperty(n,e,{value:t,enumerable:!0,configurable:!0,writable:!0}):n[e]=t,n}const Xe=n=>sr(n)?void 0:n,Wr=n=>n!==void 0;function sr(n){return n===void 0||(n instanceof Map?n.size===0:Object.values(n).filter(Wr).length===0)}function M(n){let e=new Error(n);throw delete e.stack,e}function _e(n){return(n=function(e){for(;e.endsWith("\0");)e=e.slice(0,-1);return e}(n).trim())===""?void 0:n}function gt(n){let e=function(t){let r=0;return t.ifd0.enabled&&(r+=1024),t.exif.enabled&&(r+=2048),t.makerNote&&(r+=2048),t.userComment&&(r+=1024),t.gps.enabled&&(r+=512),t.interop.enabled&&(r+=100),t.ifd1.enabled&&(r+=1024),r+2048}(n);return n.jfif.enabled&&(e+=50),n.xmp.enabled&&(e+=2e4),n.iptc.enabled&&(e+=14e3),n.icc.enabled&&(e+=6e3),e}const mt=n=>String.fromCharCode.apply(null,n),Mt=typeof TextDecoder<"u"?new TextDecoder("utf-8"):void 0;function ar(n){return Mt?Mt.decode(n):Qe?Buffer.from(n).toString("utf8"):decodeURIComponent(escape(mt(n)))}class j{static from(e,t){return e instanceof this&&e.le===t?e:new j(e,void 0,void 0,t)}constructor(e,t=0,r,i){if(typeof i=="boolean"&&(this.le=i),Array.isArray(e)&&(e=new Uint8Array(e)),e===0)this.byteOffset=0,this.byteLength=0;else if(e instanceof ArrayBuffer){r===void 0&&(r=e.byteLength-t);let s=new DataView(e,t,r);this._swapDataView(s)}else if(e instanceof Uint8Array||e instanceof DataView||e instanceof j){r===void 0&&(r=e.byteLength-t),(t+=e.byteOffset)+r>e.byteOffset+e.byteLength&&M("Creating view outside of available memory in ArrayBuffer");let s=new DataView(e.buffer,t,r);this._swapDataView(s)}else if(typeof e=="number"){let s=new DataView(new ArrayBuffer(e));this._swapDataView(s)}else M("Invalid input argument for BufferView: "+e)}_swapArrayBuffer(e){this._swapDataView(new DataView(e))}_swapBuffer(e){this._swapDataView(new DataView(e.buffer,e.byteOffset,e.byteLength))}_swapDataView(e){this.dataView=e,this.buffer=e.buffer,this.byteOffset=e.byteOffset,this.byteLength=e.byteLength}_lengthToEnd(e){return this.byteLength-e}set(e,t,r=j){return e instanceof DataView||e instanceof j?e=new Uint8Array(e.buffer,e.byteOffset,e.byteLength):e instanceof ArrayBuffer&&(e=new Uint8Array(e)),e instanceof Uint8Array||M("BufferView.set(): Invalid data argument."),this.toUint8().set(e,t),new r(this,t,e.byteLength)}subarray(e,t){return t=t||this._lengthToEnd(e),new j(this,e,t)}toUint8(){return new Uint8Array(this.buffer,this.byteOffset,this.byteLength)}getUint8Array(e,t){return new Uint8Array(this.buffer,this.byteOffset+e,t)}getString(e=0,t=this.byteLength){return ar(this.getUint8Array(e,t))}getLatin1String(e=0,t=this.byteLength){let r=this.getUint8Array(e,t);return mt(r)}getUnicodeString(e=0,t=this.byteLength){const r=[];for(let i=0;i<t&&e+i<this.byteLength;i+=2)r.push(this.getUint16(e+i));return mt(r)}getInt8(e){return this.dataView.getInt8(e)}getUint8(e){return this.dataView.getUint8(e)}getInt16(e,t=this.le){return this.dataView.getInt16(e,t)}getInt32(e,t=this.le){return this.dataView.getInt32(e,t)}getUint16(e,t=this.le){return this.dataView.getUint16(e,t)}getUint32(e,t=this.le){return this.dataView.getUint32(e,t)}getFloat32(e,t=this.le){return this.dataView.getFloat32(e,t)}getFloat64(e,t=this.le){return this.dataView.getFloat64(e,t)}getFloat(e,t=this.le){return this.dataView.getFloat32(e,t)}getDouble(e,t=this.le){return this.dataView.getFloat64(e,t)}getUintBytes(e,t,r){switch(t){case 1:return this.getUint8(e,r);case 2:return this.getUint16(e,r);case 4:return this.getUint32(e,r);case 8:return this.getUint64&&this.getUint64(e,r)}}getUint(e,t,r){switch(t){case 8:return this.getUint8(e,r);case 16:return this.getUint16(e,r);case 32:return this.getUint32(e,r);case 64:return this.getUint64&&this.getUint64(e,r)}}toString(e){return this.dataView.toString(e,this.constructor.name)}ensureChunk(){}}function _t(n,e){M(`${n} '${e}' was not loaded, try using full build of exifr.`)}class St extends Map{constructor(e){super(),this.kind=e}get(e,t){return this.has(e)||_t(this.kind,e),t&&(e in t||function(r,i){M(`Unknown ${r} '${i}'.`)}(this.kind,e),t[e].enabled||_t(this.kind,e)),super.get(e)}keyList(){return Array.from(this.keys())}}var ee=new St("file parser"),F=new St("segment parser"),te=new St("file reader");function Xr(n,e){return typeof n=="string"?Et(n,e):Fe&&!Gr&&n instanceof HTMLImageElement?Et(n.src,e):n instanceof Uint8Array||n instanceof ArrayBuffer||n instanceof DataView?new j(n):Fe&&n instanceof Blob?yt(n,e,"blob",xe):void M("Invalid input argument")}function Et(n,e){return(t=n).startsWith("data:")||t.length>1e4?wt(n,e,"base64"):je&&n.includes("://")?yt(n,e,"url",be):je?wt(n,e,"fs"):Fe?yt(n,e,"url",be):void M("Invalid input argument");var t}async function yt(n,e,t,r){return te.has(t)?wt(n,e,t):r?async function(i,s){let a=await s(i);return new j(a)}(n,r):void M(`Parser ${t} is not loaded`)}async function wt(n,e,t){let r=new(te.get(t))(n,e);return await r.read(),r}const be=n=>xt(n).then(e=>e.arrayBuffer()),xe=n=>new Promise((e,t)=>{let r=new FileReader;r.onloadend=()=>e(r.result||new ArrayBuffer),r.onerror=t,r.readAsArrayBuffer(n)});class Kr extends Map{get tagKeys(){return this.allKeys||(this.allKeys=Array.from(this.keys())),this.allKeys}get tagValues(){return this.allValues||(this.allValues=Array.from(this.values())),this.allValues}}function L(n,e,t){let r=new Kr;for(let[i,s]of t)r.set(i,s);if(Array.isArray(e))for(let i of e)n.set(i,r);else n.set(e,r);return r}function Se(n,e,t){let r,i=n.get(e);for(r of t)i.set(r[0],r[1])}const N=new Map,$=new Map,he=new Map,oe=["chunked","firstChunkSize","firstChunkSizeNode","firstChunkSizeBrowser","chunkSize","chunkLimit"],Ae=["jfif","xmp","icc","iptc","ihdr"],ke=["tiff",...Ae],I=["ifd0","ifd1","exif","gps","interop"],le=[...ke,...I],de=["makerNote","userComment"],Ie=["translateKeys","translateValues","reviveValues","multiSegment"],ce=[...Ie,"sanitize","mergeOutput","silentErrors"];class or{get translate(){return this.translateKeys||this.translateValues||this.reviveValues}}class Le extends or{get needed(){return this.enabled||this.deps.size>0}constructor(e,t,r,i){if(super(),m(this,"enabled",!1),m(this,"skip",new Set),m(this,"pick",new Set),m(this,"deps",new Set),m(this,"translateKeys",!1),m(this,"translateValues",!1),m(this,"reviveValues",!1),this.key=e,this.enabled=t,this.parse=this.enabled,this.applyInheritables(i),this.canBeFiltered=I.includes(e),this.canBeFiltered&&(this.dict=N.get(e)),r!==void 0)if(Array.isArray(r))this.parse=this.enabled=!0,this.canBeFiltered&&r.length>0&&this.translateTagSet(r,this.pick);else if(typeof r=="object"){if(this.enabled=!0,this.parse=r.parse!==!1,this.canBeFiltered){let{pick:s,skip:a}=r;s&&s.length>0&&this.translateTagSet(s,this.pick),a&&a.length>0&&this.translateTagSet(a,this.skip)}this.applyInheritables(r)}else r===!0||r===!1?this.parse=this.enabled=r:M(`Invalid options argument: ${r}`)}applyInheritables(e){let t,r;for(t of Ie)r=e[t],r!==void 0&&(this[t]=r)}translateTagSet(e,t){if(this.dict){let r,i,{tagKeys:s,tagValues:a}=this.dict;for(r of e)typeof r=="string"?(i=a.indexOf(r),i===-1&&(i=s.indexOf(Number(r))),i!==-1&&t.add(Number(s[i]))):t.add(r)}else for(let r of e)t.add(r)}finalizeFilters(){!this.enabled&&this.deps.size>0?(this.enabled=!0,Ke(this.pick,this.deps)):this.enabled&&this.pick.size>0&&Ke(this.pick,this.deps)}}var G={jfif:!1,tiff:!0,xmp:!1,icc:!1,iptc:!1,ifd0:!0,ifd1:!1,exif:!0,gps:!0,interop:!1,ihdr:void 0,makerNote:!1,userComment:!1,multiSegment:!1,skip:[],pick:[],translateKeys:!0,translateValues:!0,reviveValues:!0,sanitize:!0,mergeOutput:!0,silentErrors:!0,chunked:!0,firstChunkSize:void 0,firstChunkSizeNode:512,firstChunkSizeBrowser:65536,chunkSize:65536,chunkLimit:5},Ut=new Map;class Ce extends or{static useCached(e){let t=Ut.get(e);return t!==void 0||(t=new this(e),Ut.set(e,t)),t}constructor(e){super(),e===!0?this.setupFromTrue():e===void 0?this.setupFromUndefined():Array.isArray(e)?this.setupFromArray(e):typeof e=="object"?this.setupFromObject(e):M(`Invalid options argument ${e}`),this.firstChunkSize===void 0&&(this.firstChunkSize=Fe?this.firstChunkSizeBrowser:this.firstChunkSizeNode),this.mergeOutput&&(this.ifd1.enabled=!1),this.filterNestedSegmentTags(),this.traverseTiffDependencyTree(),this.checkLoadedPlugins()}setupFromUndefined(){let e;for(e of oe)this[e]=G[e];for(e of ce)this[e]=G[e];for(e of de)this[e]=G[e];for(e of le)this[e]=new Le(e,G[e],void 0,this)}setupFromTrue(){let e;for(e of oe)this[e]=G[e];for(e of ce)this[e]=G[e];for(e of de)this[e]=!0;for(e of le)this[e]=new Le(e,!0,void 0,this)}setupFromArray(e){let t;for(t of oe)this[t]=G[t];for(t of ce)this[t]=G[t];for(t of de)this[t]=G[t];for(t of le)this[t]=new Le(t,!1,void 0,this);this.setupGlobalFilters(e,void 0,I)}setupFromObject(e){let t;for(t of(I.ifd0=I.ifd0||I.image,I.ifd1=I.ifd1||I.thumbnail,Object.assign(this,e),oe))this[t]=dt(e[t],G[t]);for(t of ce)this[t]=dt(e[t],G[t]);for(t of de)this[t]=dt(e[t],G[t]);for(t of ke)this[t]=new Le(t,G[t],e[t],this);for(t of I)this[t]=new Le(t,G[t],e[t],this.tiff);this.setupGlobalFilters(e.pick,e.skip,I,le),e.tiff===!0?this.batchEnableWithBool(I,!0):e.tiff===!1?this.batchEnableWithUserValue(I,e):Array.isArray(e.tiff)?this.setupGlobalFilters(e.tiff,void 0,I):typeof e.tiff=="object"&&this.setupGlobalFilters(e.tiff.pick,e.tiff.skip,I)}batchEnableWithBool(e,t){for(let r of e)this[r].enabled=t}batchEnableWithUserValue(e,t){for(let r of e){let i=t[r];this[r].enabled=i!==!1&&i!==void 0}}setupGlobalFilters(e,t,r,i=r){if(e&&e.length){for(let a of i)this[a].enabled=!1;let s=Nt(e,r);for(let[a,o]of s)Ke(this[a].pick,o),this[a].enabled=!0}else if(t&&t.length){let s=Nt(t,r);for(let[a,o]of s)Ke(this[a].skip,o)}}filterNestedSegmentTags(){let{ifd0:e,exif:t,xmp:r,iptc:i,icc:s}=this;this.makerNote?t.deps.add(37500):t.skip.add(37500),this.userComment?t.deps.add(37510):t.skip.add(37510),r.enabled||e.skip.add(700),i.enabled||e.skip.add(33723),s.enabled||e.skip.add(34675)}traverseTiffDependencyTree(){let{ifd0:e,exif:t,gps:r,interop:i}=this;i.needed&&(t.deps.add(40965),e.deps.add(40965)),t.needed&&e.deps.add(34665),r.needed&&e.deps.add(34853),this.tiff.enabled=I.some(s=>this[s].enabled===!0)||this.makerNote||this.userComment;for(let s of I)this[s].finalizeFilters()}get onlyTiff(){return!Ae.map(e=>this[e].enabled).some(e=>e===!0)&&this.tiff.enabled}checkLoadedPlugins(){for(let e of ke)this[e].enabled&&!F.has(e)&&_t("segment parser",e)}}function Nt(n,e){let t,r,i,s,a=[];for(i of e){for(s of(t=N.get(i),r=[],t))(n.includes(s[0])||n.includes(s[1]))&&r.push(s[0]);r.length&&a.push([i,r])}return a}function dt(n,e){return n!==void 0?n:e!==void 0?e:void 0}function Ke(n,e){for(let t of e)n.add(t)}m(Ce,"default",G);class pe{constructor(e){m(this,"parsers",{}),m(this,"output",{}),m(this,"errors",[]),m(this,"pushToErrors",t=>this.errors.push(t)),this.options=Ce.useCached(e)}async read(e){this.file=await Xr(e,this.options)}setup(){if(this.fileParser)return;let{file:e}=this,t=e.getUint16(0);for(let[r,i]of ee)if(i.canHandle(e,t))return this.fileParser=new i(this.options,this.file,this.parsers),e[r]=!0;this.file.close&&this.file.close(),M("Unknown file format")}async parse(){let{output:e,errors:t}=this;return this.setup(),this.options.silentErrors?(await this.executeParsers().catch(this.pushToErrors),t.push(...this.fileParser.errors)):await this.executeParsers(),this.file.close&&this.file.close(),this.options.silentErrors&&t.length>0&&(e.errors=t),Xe(e)}async executeParsers(){let{output:e}=this;await this.fileParser.parse();let t=Object.values(this.parsers).map(async r=>{let i=await r.parse();r.assignToOutput(e,i)});this.options.silentErrors&&(t=t.map(r=>r.catch(this.pushToErrors))),await Promise.all(t)}async extractThumbnail(){this.setup();let{options:e,file:t}=this,r=F.get("tiff",e);var i;if(t.tiff?i={start:0,type:"tiff"}:t.jpeg&&(i=await this.fileParser.getOrFindSegment("tiff")),i===void 0)return;let s=await this.fileParser.ensureSegmentChunk(i),a=this.parsers.tiff=new r(s,e,t),o=await a.extractThumbnail();return t.close&&t.close(),o}}async function Ze(n,e){let t=new pe(e);return await t.read(n),t.parse()}var qr=Object.freeze({__proto__:null,parse:Ze,Exifr:pe,fileParsers:ee,segmentParsers:F,fileReaders:te,tagKeys:N,tagValues:$,tagRevivers:he,createDictionary:L,extendDictionary:Se,fetchUrlAsArrayBuffer:be,readBlobAsArrayBuffer:xe,chunkedProps:oe,otherSegments:Ae,segments:ke,tiffBlocks:I,segmentsAndBlocks:le,tiffExtractables:de,inheritables:Ie,allFormatters:ce,Options:Ce});class et{constructor(e,t,r){m(this,"errors",[]),m(this,"ensureSegmentChunk",async i=>{let s=i.start,a=i.size||65536;if(this.file.chunked)if(this.file.available(s,a))i.chunk=this.file.subarray(s,a);else try{i.chunk=await this.file.readChunk(s,a)}catch(o){M(`Couldn't read segment: ${JSON.stringify(i)}. ${o.message}`)}else this.file.byteLength>s+a?i.chunk=this.file.subarray(s,a):i.size===void 0?i.chunk=this.file.subarray(s):M("Segment unreachable: "+JSON.stringify(i));return i.chunk}),this.extendOptions&&this.extendOptions(e),this.options=e,this.file=t,this.parsers=r}injectSegment(e,t){this.options[e].enabled&&this.createParser(e,t)}createParser(e,t){let r=new(F.get(e))(t,this.options,this.file);return this.parsers[e]=r}createParsers(e){for(let t of e){let{type:r,chunk:i}=t,s=this.options[r];if(s&&s.enabled){let a=this.parsers[r];a&&a.append||a||this.createParser(r,i)}}}async readSegments(e){let t=e.map(this.ensureSegmentChunk);await Promise.all(t)}}class Y{static findPosition(e,t){let r=e.getUint16(t+2)+2,i=typeof this.headerLength=="function"?this.headerLength(e,t,r):this.headerLength,s=t+i,a=r-i;return{offset:t,length:r,headerLength:i,start:s,size:a,end:s+a}}static parse(e,t={}){return new this(e,new Ce({[this.type]:t}),e).parse()}normalizeInput(e){return e instanceof j?e:new j(e)}constructor(e,t={},r){m(this,"errors",[]),m(this,"raw",new Map),m(this,"handleError",i=>{if(!this.options.silentErrors)throw i;this.errors.push(i.message)}),this.chunk=this.normalizeInput(e),this.file=r,this.type=this.constructor.type,this.globalOptions=this.options=t,this.localOptions=t[this.type],this.canTranslate=this.localOptions&&this.localOptions.translate}translate(){this.canTranslate&&(this.translated=this.translateBlock(this.raw,this.type))}get output(){return this.translated?this.translated:this.raw?Object.fromEntries(this.raw):void 0}translateBlock(e,t){let r=he.get(t),i=$.get(t),s=N.get(t),a=this.options[t],o=a.reviveValues&&!!r,l=a.translateValues&&!!i,h=a.translateKeys&&!!s,c={};for(let[d,u]of e)o&&r.has(d)?u=r.get(d)(u):l&&i.has(d)&&(u=this.translateValue(u,i.get(d))),h&&s.has(d)&&(d=s.get(d)||d),c[d]=u;return c}translateValue(e,t){return t[e]||t.DEFAULT||e}assignToOutput(e,t){this.assignObjectToOutput(e,this.constructor.type,t)}assignObjectToOutput(e,t,r){if(this.globalOptions.mergeOutput)return Object.assign(e,r);e[t]?Object.assign(e[t],r):e[t]=r}}m(Y,"headerLength",4),m(Y,"type",void 0),m(Y,"multiSegment",!1),m(Y,"canHandle",()=>!1);function Yr(n){return n===192||n===194||n===196||n===219||n===221||n===218||n===254}function $r(n){return n>=224&&n<=239}function Jr(n,e,t){for(let[r,i]of F)if(i.canHandle(n,e,t))return r}class Bt extends et{constructor(...e){super(...e),m(this,"appSegments",[]),m(this,"jpegSegments",[]),m(this,"unknownSegments",[])}static canHandle(e,t){return t===65496}async parse(){await this.findAppSegments(),await this.readSegments(this.appSegments),this.mergeMultiSegments(),this.createParsers(this.mergedAppSegments||this.appSegments)}setupSegmentFinderArgs(e){e===!0?(this.findAll=!0,this.wanted=new Set(F.keyList())):(e=e===void 0?F.keyList().filter(t=>this.options[t].enabled):e.filter(t=>this.options[t].enabled&&F.has(t)),this.findAll=!1,this.remaining=new Set(e),this.wanted=new Set(e)),this.unfinishedMultiSegment=!1}async findAppSegments(e=0,t){this.setupSegmentFinderArgs(t);let{file:r,findAll:i,wanted:s,remaining:a}=this;if(!i&&this.file.chunked&&(i=Array.from(s).some(o=>{let l=F.get(o),h=this.options[o];return l.multiSegment&&h.multiSegment}),i&&await this.file.readWhole()),e=this.findAppSegmentsInRange(e,r.byteLength),!this.options.onlyTiff&&r.chunked){let o=!1;for(;a.size>0&&!o&&(r.canReadNextChunk||this.unfinishedMultiSegment);){let{nextChunkOffset:l}=r,h=this.appSegments.some(c=>!this.file.available(c.offset||c.start,c.length||c.size));if(o=e>l&&!h?!await r.readNextChunk(e):!await r.readNextChunk(l),(e=this.findAppSegmentsInRange(e,r.byteLength))===void 0)return}}}findAppSegmentsInRange(e,t){t-=2;let r,i,s,a,o,l,{file:h,findAll:c,wanted:d,remaining:u,options:p}=this;for(;e<t;e++)if(h.getUint8(e)===255){if(r=h.getUint8(e+1),$r(r)){if(i=h.getUint16(e+2),s=Jr(h,e,i),s&&d.has(s)&&(a=F.get(s),o=a.findPosition(h,e),l=p[s],o.type=s,this.appSegments.push(o),!c&&(a.multiSegment&&l.multiSegment?(this.unfinishedMultiSegment=o.chunkNumber<o.chunkCount,this.unfinishedMultiSegment||u.delete(s)):u.delete(s),u.size===0)))break;p.recordUnknownSegments&&(o=Y.findPosition(h,e),o.marker=r,this.unknownSegments.push(o)),e+=i+1}else if(Yr(r)){if(i=h.getUint16(e+2),r===218&&p.stopAfterSos!==!1)return;p.recordJpegSegments&&this.jpegSegments.push({offset:e,length:i,marker:r}),e+=i+1}}return e}mergeMultiSegments(){if(!this.appSegments.some(t=>t.multiSegment))return;let e=function(t,r){let i,s,a,o=new Map;for(let l=0;l<t.length;l++)i=t[l],s=i[r],o.has(s)?a=o.get(s):o.set(s,a=[]),a.push(i);return Array.from(o)}(this.appSegments,"type");this.mergedAppSegments=e.map(([t,r])=>{let i=F.get(t,this.options);return i.handleMultiSegments?{type:t,chunk:i.handleMultiSegments(r)}:r[0]})}getSegment(e){return this.appSegments.find(t=>t.type===e)}async getOrFindSegment(e){let t=this.getSegment(e);return t===void 0&&(await this.findAppSegments(0,[e]),t=this.getSegment(e)),t}}m(Bt,"type","jpeg"),ee.set("jpeg",Bt);const Qr=[void 0,1,1,2,4,8,1,1,2,4,8,4,8,4];class Zr extends Y{parseHeader(){var e=this.chunk.getUint16();e===18761?this.le=!0:e===19789&&(this.le=!1),this.chunk.le=this.le,this.headerParsed=!0}parseTags(e,t,r=new Map){let{pick:i,skip:s}=this.options[t];i=new Set(i);let a=i.size>0,o=s.size===0,l=this.chunk.getUint16(e);e+=2;for(let h=0;h<l;h++){let c=this.chunk.getUint16(e);if(a){if(i.has(c)&&(r.set(c,this.parseTag(e,c,t)),i.delete(c),i.size===0))break}else!o&&s.has(c)||r.set(c,this.parseTag(e,c,t));e+=12}return r}parseTag(e,t,r){let{chunk:i}=this,s=i.getUint16(e+2),a=i.getUint32(e+4),o=Qr[s];if(o*a<=4?e+=8:e=i.getUint32(e+8),(s<1||s>13)&&M(`Invalid TIFF value type. block: ${r.toUpperCase()}, tag: ${t.toString(16)}, type: ${s}, offset ${e}`),e>i.byteLength&&M(`Invalid TIFF value offset. block: ${r.toUpperCase()}, tag: ${t.toString(16)}, type: ${s}, offset ${e} is outside of chunk size ${i.byteLength}`),s===1)return i.getUint8Array(e,a);if(s===2)return _e(i.getString(e,a));if(s===7)return i.getUint8Array(e,a);if(a===1)return this.parseTagValue(s,e);{let l=new(function(c){switch(c){case 1:return Uint8Array;case 3:return Uint16Array;case 4:return Uint32Array;case 5:return Array;case 6:return Int8Array;case 8:return Int16Array;case 9:return Int32Array;case 10:return Array;case 11:return Float32Array;case 12:return Float64Array;default:return Array}}(s))(a),h=o;for(let c=0;c<a;c++)l[c]=this.parseTagValue(s,e),e+=h;return l}}parseTagValue(e,t){let{chunk:r}=this;switch(e){case 1:return r.getUint8(t);case 3:return r.getUint16(t);case 4:return r.getUint32(t);case 5:return r.getUint32(t)/r.getUint32(t+4);case 6:return r.getInt8(t);case 8:return r.getInt16(t);case 9:return r.getInt32(t);case 10:return r.getInt32(t)/r.getInt32(t+4);case 11:return r.getFloat(t);case 12:return r.getDouble(t);case 13:return r.getUint32(t);default:M(`Invalid tiff type ${e}`)}}}class ct extends Zr{static canHandle(e,t){return e.getUint8(t+1)===225&&e.getUint32(t+4)===1165519206&&e.getUint16(t+8)===0}async parse(){this.parseHeader();let{options:e}=this;return e.ifd0.enabled&&await this.parseIfd0Block(),e.exif.enabled&&await this.safeParse("parseExifBlock"),e.gps.enabled&&await this.safeParse("parseGpsBlock"),e.interop.enabled&&await this.safeParse("parseInteropBlock"),e.ifd1.enabled&&await this.safeParse("parseThumbnailBlock"),this.createOutput()}safeParse(e){let t=this[e]();return t.catch!==void 0&&(t=t.catch(this.handleError)),t}findIfd0Offset(){this.ifd0Offset===void 0&&(this.ifd0Offset=this.chunk.getUint32(4))}findIfd1Offset(){if(this.ifd1Offset===void 0){this.findIfd0Offset();let e=this.chunk.getUint16(this.ifd0Offset),t=this.ifd0Offset+2+12*e;this.ifd1Offset=this.chunk.getUint32(t)}}parseBlock(e,t){let r=new Map;return this[t]=r,this.parseTags(e,t,r),r}async parseIfd0Block(){if(this.ifd0)return;let{file:e}=this;this.findIfd0Offset(),this.ifd0Offset<8&&M("Malformed EXIF data"),!e.chunked&&this.ifd0Offset>e.byteLength&&M(`IFD0 offset points to outside of file.
this.ifd0Offset: ${this.ifd0Offset}, file.byteLength: ${e.byteLength}`),e.tiff&&await e.ensureChunk(this.ifd0Offset,gt(this.options));let t=this.parseBlock(this.ifd0Offset,"ifd0");return t.size!==0?(this.exifOffset=t.get(34665),this.interopOffset=t.get(40965),this.gpsOffset=t.get(34853),this.xmp=t.get(700),this.iptc=t.get(33723),this.icc=t.get(34675),this.options.sanitize&&(t.delete(34665),t.delete(40965),t.delete(34853),t.delete(700),t.delete(33723),t.delete(34675)),t):void 0}async parseExifBlock(){if(this.exif||(this.ifd0||await this.parseIfd0Block(),this.exifOffset===void 0))return;this.file.tiff&&await this.file.ensureChunk(this.exifOffset,gt(this.options));let e=this.parseBlock(this.exifOffset,"exif");return this.interopOffset||(this.interopOffset=e.get(40965)),this.makerNote=e.get(37500),this.userComment=e.get(37510),this.options.sanitize&&(e.delete(40965),e.delete(37500),e.delete(37510)),this.unpack(e,41728),this.unpack(e,41729),e}unpack(e,t){let r=e.get(t);r&&r.length===1&&e.set(t,r[0])}async parseGpsBlock(){if(this.gps||(this.ifd0||await this.parseIfd0Block(),this.gpsOffset===void 0))return;let e=this.parseBlock(this.gpsOffset,"gps");return e&&e.has(2)&&e.has(4)&&(e.set("latitude",Vt(...e.get(2),e.get(1))),e.set("longitude",Vt(...e.get(4),e.get(3)))),e}async parseInteropBlock(){if(!this.interop&&(this.ifd0||await this.parseIfd0Block(),this.interopOffset!==void 0||this.exif||await this.parseExifBlock(),this.interopOffset!==void 0))return this.parseBlock(this.interopOffset,"interop")}async parseThumbnailBlock(e=!1){if(!this.ifd1&&!this.ifd1Parsed&&(!this.options.mergeOutput||e))return this.findIfd1Offset(),this.ifd1Offset>0&&(this.parseBlock(this.ifd1Offset,"ifd1"),this.ifd1Parsed=!0),this.ifd1}async extractThumbnail(){if(this.headerParsed||this.parseHeader(),this.ifd1Parsed||await this.parseThumbnailBlock(!0),this.ifd1===void 0)return;let e=this.ifd1.get(513),t=this.ifd1.get(514);return this.chunk.getUint8Array(e,t)}get image(){return this.ifd0}get thumbnail(){return this.ifd1}createOutput(){let e,t,r,i={};for(t of I)if(e=this[t],!sr(e))if(r=this.canTranslate?this.translateBlock(e,t):Object.fromEntries(e),this.options.mergeOutput){if(t==="ifd1")continue;Object.assign(i,r)}else i[t]=r;return this.makerNote&&(i.makerNote=this.makerNote),this.userComment&&(i.userComment=this.userComment),i}assignToOutput(e,t){if(this.globalOptions.mergeOutput)Object.assign(e,t);else for(let[r,i]of Object.entries(t))this.assignObjectToOutput(e,r,i)}}function Vt(n,e,t,r){var i=n+e/60+t/3600;return r!=="S"&&r!=="W"||(i*=-1),i}m(ct,"type","tiff"),m(ct,"headerLength",10),F.set("tiff",ct);var en=Object.freeze({__proto__:null,default:qr,Exifr:pe,fileParsers:ee,segmentParsers:F,fileReaders:te,tagKeys:N,tagValues:$,tagRevivers:he,createDictionary:L,extendDictionary:Se,fetchUrlAsArrayBuffer:be,readBlobAsArrayBuffer:xe,chunkedProps:oe,otherSegments:Ae,segments:ke,tiffBlocks:I,segmentsAndBlocks:le,tiffExtractables:de,inheritables:Ie,allFormatters:ce,Options:Ce,parse:Ze});const kt={ifd0:!1,ifd1:!1,exif:!1,gps:!1,interop:!1,sanitize:!1,reviveValues:!0,translateKeys:!1,translateValues:!1,mergeOutput:!1},Ct=Object.assign({},kt,{firstChunkSize:4e4,gps:[1,2,3,4]});async function lr(n){let e=new pe(Ct);await e.read(n);let t=await e.parse();if(t&&t.gps){let{latitude:r,longitude:i}=t.gps;return{latitude:r,longitude:i}}}const vt=Object.assign({},kt,{tiff:!1,ifd1:!0,mergeOutput:!1});async function dr(n){let e=new pe(vt);await e.read(n);let t=await e.extractThumbnail();return t&&Qe?Je.from(t):t}async function cr(n){let e=await this.thumbnail(n);if(e!==void 0){let t=new Blob([e]);return URL.createObjectURL(t)}}const Pt=Object.assign({},kt,{firstChunkSize:4e4,ifd0:[274]});async function At(n){let e=new pe(Pt);await e.read(n);let t=await e.parse();if(t&&t.ifd0)return t.ifd0[274]}const It=Object.freeze({1:{dimensionSwapped:!1,scaleX:1,scaleY:1,deg:0,rad:0},2:{dimensionSwapped:!1,scaleX:-1,scaleY:1,deg:0,rad:0},3:{dimensionSwapped:!1,scaleX:1,scaleY:1,deg:180,rad:180*Math.PI/180},4:{dimensionSwapped:!1,scaleX:-1,scaleY:1,deg:180,rad:180*Math.PI/180},5:{dimensionSwapped:!0,scaleX:1,scaleY:-1,deg:90,rad:90*Math.PI/180},6:{dimensionSwapped:!0,scaleX:1,scaleY:1,deg:90,rad:90*Math.PI/180},7:{dimensionSwapped:!0,scaleX:1,scaleY:-1,deg:270,rad:270*Math.PI/180},8:{dimensionSwapped:!0,scaleX:1,scaleY:1,deg:270,rad:270*Math.PI/180}});let ye=!0,we=!0;if(typeof navigator=="object"){let n=navigator.userAgent;if(n.includes("iPad")||n.includes("iPhone")){let e=n.match(/OS (\d+)_(\d+)/);if(e){let[,t,r]=e;ye=Number(t)+.1*Number(r)<13.4,we=!1}}else if(n.includes("OS X 10")){let[,e]=n.match(/OS X 10[_.](\d+)/);ye=we=Number(e)<15}if(n.includes("Chrome/")){let[,e]=n.match(/Chrome\/(\d+)/);ye=we=Number(e)<81}else if(n.includes("Firefox/")){let[,e]=n.match(/Firefox\/(\d+)/);ye=we=Number(e)<77}}async function hr(n){let e=await At(n);return Object.assign({canvas:ye,css:we},It[e])}class tn extends j{constructor(...e){super(...e),m(this,"ranges",new rn),this.byteLength!==0&&this.ranges.add(0,this.byteLength)}_tryExtend(e,t,r){if(e===0&&this.byteLength===0&&r){let i=new DataView(r.buffer||r,r.byteOffset,r.byteLength);this._swapDataView(i)}else{let i=e+t;if(i>this.byteLength){let{dataView:s}=this._extend(i);this._swapDataView(s)}}}_extend(e){let t;t=Qe?Je.allocUnsafe(e):new Uint8Array(e);let r=new DataView(t.buffer,t.byteOffset,t.byteLength);return t.set(new Uint8Array(this.buffer,this.byteOffset,this.byteLength),0),{uintView:t,dataView:r}}subarray(e,t,r=!1){return t=t||this._lengthToEnd(e),r&&this._tryExtend(e,t),this.ranges.add(e,t),super.subarray(e,t)}set(e,t,r=!1){r&&this._tryExtend(t,e.byteLength,e);let i=super.set(e,t);return this.ranges.add(t,i.byteLength),i}async ensureChunk(e,t){this.chunked&&(this.ranges.available(e,t)||await this.readChunk(e,t))}available(e,t){return this.ranges.available(e,t)}}class rn{constructor(){m(this,"list",[])}get length(){return this.list.length}add(e,t,r=0){let i=e+t,s=this.list.filter(a=>zt(e,a.offset,i)||zt(e,a.end,i));if(s.length>0){e=Math.min(e,...s.map(o=>o.offset)),i=Math.max(i,...s.map(o=>o.end)),t=i-e;let a=s.shift();a.offset=e,a.length=t,a.end=i,this.list=this.list.filter(o=>!s.includes(o))}else this.list.push({offset:e,length:t,end:i})}available(e,t){let r=e+t;return this.list.some(i=>i.offset<=e&&r<=i.end)}}function zt(n,e,t){return n<=e&&e<=t}class tt extends tn{constructor(e,t){super(0),m(this,"chunksRead",0),this.input=e,this.options=t}async readWhole(){this.chunked=!1,await this.readChunk(this.nextChunkOffset)}async readChunked(){this.chunked=!0,await this.readChunk(0,this.options.firstChunkSize)}async readNextChunk(e=this.nextChunkOffset){if(this.fullyRead)return this.chunksRead++,!1;let t=this.options.chunkSize,r=await this.readChunk(e,t);return!!r&&r.byteLength===t}async readChunk(e,t){if(this.chunksRead++,(t=this.safeWrapAddress(e,t))!==0)return this._readChunk(e,t)}safeWrapAddress(e,t){return this.size!==void 0&&e+t>this.size?Math.max(0,this.size-e):t}get nextChunkOffset(){if(this.ranges.list.length!==0)return this.ranges.list[0].length}get canReadNextChunk(){return this.chunksRead<this.options.chunkLimit}get fullyRead(){return this.size!==void 0&&this.nextChunkOffset===this.size}read(){return this.options.chunked?this.readChunked():this.readWhole()}close(){}}te.set("blob",class extends tt{async readWhole(){this.chunked=!1;let n=await xe(this.input);this._swapArrayBuffer(n)}readChunked(){return this.chunked=!0,this.size=this.input.size,super.readChunked()}async _readChunk(n,e){let t=e?n+e:void 0,r=this.input.slice(n,t),i=await xe(r);return this.set(i,n,!0)}});var nn=Object.freeze({__proto__:null,default:en,Exifr:pe,fileParsers:ee,segmentParsers:F,fileReaders:te,tagKeys:N,tagValues:$,tagRevivers:he,createDictionary:L,extendDictionary:Se,fetchUrlAsArrayBuffer:be,readBlobAsArrayBuffer:xe,chunkedProps:oe,otherSegments:Ae,segments:ke,tiffBlocks:I,segmentsAndBlocks:le,tiffExtractables:de,inheritables:Ie,allFormatters:ce,Options:Ce,parse:Ze,gpsOnlyOptions:Ct,gps:lr,thumbnailOnlyOptions:vt,thumbnail:dr,thumbnailUrl:cr,orientationOnlyOptions:Pt,orientation:At,rotations:It,get rotateCanvas(){return ye},get rotateCss(){return we},rotation:hr});te.set("url",class extends tt{async readWhole(){this.chunked=!1;let n=await be(this.input);n instanceof ArrayBuffer?this._swapArrayBuffer(n):n instanceof Uint8Array&&this._swapBuffer(n)}async _readChunk(n,e){let t=e?n+e-1:void 0,r=this.options.httpHeaders||{};(n||t)&&(r.range=`bytes=${[n,t].join("-")}`);let i=await xt(this.input,{headers:r}),s=await i.arrayBuffer(),a=s.byteLength;if(i.status!==416)return a!==e&&(this.size=n+a),this.set(s,n,!0)}});j.prototype.getUint64=function(n){let e=this.getUint32(n),t=this.getUint32(n+4);return e<1048575?e<<32|t:typeof Ne!==void 0?(console.warn("Using BigInt because of type 64uint but JS can only handle 53b numbers."),Ne(e)<<Ne(32)|Ne(t)):void M("Trying to read 64b value but JS can only handle 53b numbers.")};class sn extends et{parseBoxes(e=0){let t=[];for(;e<this.file.byteLength-4;){let r=this.parseBoxHead(e);if(t.push(r),r.length===0)break;e+=r.length}return t}parseSubBoxes(e){e.boxes=this.parseBoxes(e.start)}findBox(e,t){return e.boxes===void 0&&this.parseSubBoxes(e),e.boxes.find(r=>r.kind===t)}parseBoxHead(e){let t=this.file.getUint32(e),r=this.file.getString(e+4,4),i=e+8;return t===1&&(t=this.file.getUint64(e+8),i+=8),{offset:e,length:t,kind:r,start:i}}parseBoxFullHead(e){if(e.version!==void 0)return;let t=this.file.getUint32(e.start);e.version=t>>24,e.start+=4}}class pr extends sn{static canHandle(e,t){if(t!==0)return!1;let r=e.getUint16(2);if(r>50)return!1;let i=16,s=[];for(;i<r;)s.push(e.getString(i,4)),i+=4;return s.includes(this.type)}async parse(){let e=this.file.getUint32(0),t=this.parseBoxHead(e);for(;t.kind!=="meta";)e+=t.length,await this.file.ensureChunk(e,16),t=this.parseBoxHead(e);await this.file.ensureChunk(t.offset,t.length),this.parseBoxFullHead(t),this.parseSubBoxes(t),this.options.icc.enabled&&await this.findIcc(t),this.options.tiff.enabled&&await this.findExif(t)}async registerSegment(e,t,r){await this.file.ensureChunk(t,r);let i=this.file.subarray(t,r);this.createParser(e,i)}async findIcc(e){let t=this.findBox(e,"iprp");if(t===void 0)return;let r=this.findBox(t,"ipco");if(r===void 0)return;let i=this.findBox(r,"colr");i!==void 0&&await this.registerSegment("icc",i.offset+12,i.length)}async findExif(e){let t=this.findBox(e,"iinf");if(t===void 0)return;let r=this.findBox(e,"iloc");if(r===void 0)return;let i=this.findExifLocIdInIinf(t),s=this.findExtentInIloc(r,i);if(s===void 0)return;let[a,o]=s;await this.file.ensureChunk(a,o);let l=4+this.file.getUint32(a);a+=l,o-=l,await this.registerSegment("tiff",a,o)}findExifLocIdInIinf(e){this.parseBoxFullHead(e);let t,r,i,s,a=e.start,o=this.file.getUint16(a);for(a+=2;o--;){if(t=this.parseBoxHead(a),this.parseBoxFullHead(t),r=t.start,t.version>=2&&(i=t.version===3?4:2,s=this.file.getString(r+i+2,4),s==="Exif"))return this.file.getUintBytes(r,i);a+=t.length}}get8bits(e){let t=this.file.getUint8(e);return[t>>4,15&t]}findExtentInIloc(e,t){this.parseBoxFullHead(e);let r=e.start,[i,s]=this.get8bits(r++),[a,o]=this.get8bits(r++),l=e.version===2?4:2,h=e.version===1||e.version===2?2:0,c=o+i+s,d=e.version===2?4:2,u=this.file.getUintBytes(r,d);for(r+=d;u--;){let p=this.file.getUintBytes(r,l);r+=l+h+2+a;let _=this.file.getUint16(r);if(r+=2,p===t)return _>1&&console.warn(`ILOC box has more than one extent but we're only processing one
Please create an issue at https://github.com/MikeKovarik/exifr with this file`),[this.file.getUintBytes(r+o,i),this.file.getUintBytes(r+o+i,s)];r+=_*c}}}class ur extends pr{}m(ur,"type","heic");class Gt extends pr{}m(Gt,"type","avif"),ee.set("heic",ur),ee.set("avif",Gt),L(N,["ifd0","ifd1"],[[256,"ImageWidth"],[257,"ImageHeight"],[258,"BitsPerSample"],[259,"Compression"],[262,"PhotometricInterpretation"],[270,"ImageDescription"],[271,"Make"],[272,"Model"],[273,"StripOffsets"],[274,"Orientation"],[277,"SamplesPerPixel"],[278,"RowsPerStrip"],[279,"StripByteCounts"],[282,"XResolution"],[283,"YResolution"],[284,"PlanarConfiguration"],[296,"ResolutionUnit"],[301,"TransferFunction"],[305,"Software"],[306,"ModifyDate"],[315,"Artist"],[316,"HostComputer"],[317,"Predictor"],[318,"WhitePoint"],[319,"PrimaryChromaticities"],[513,"ThumbnailOffset"],[514,"ThumbnailLength"],[529,"YCbCrCoefficients"],[530,"YCbCrSubSampling"],[531,"YCbCrPositioning"],[532,"ReferenceBlackWhite"],[700,"ApplicationNotes"],[33432,"Copyright"],[33723,"IPTC"],[34665,"ExifIFD"],[34675,"ICC"],[34853,"GpsIFD"],[330,"SubIFD"],[40965,"InteropIFD"],[40091,"XPTitle"],[40092,"XPComment"],[40093,"XPAuthor"],[40094,"XPKeywords"],[40095,"XPSubject"]]),L(N,"exif",[[33434,"ExposureTime"],[33437,"FNumber"],[34850,"ExposureProgram"],[34852,"SpectralSensitivity"],[34855,"ISO"],[34858,"TimeZoneOffset"],[34859,"SelfTimerMode"],[34864,"SensitivityType"],[34865,"StandardOutputSensitivity"],[34866,"RecommendedExposureIndex"],[34867,"ISOSpeed"],[34868,"ISOSpeedLatitudeyyy"],[34869,"ISOSpeedLatitudezzz"],[36864,"ExifVersion"],[36867,"DateTimeOriginal"],[36868,"CreateDate"],[36873,"GooglePlusUploadCode"],[36880,"OffsetTime"],[36881,"OffsetTimeOriginal"],[36882,"OffsetTimeDigitized"],[37121,"ComponentsConfiguration"],[37122,"CompressedBitsPerPixel"],[37377,"ShutterSpeedValue"],[37378,"ApertureValue"],[37379,"BrightnessValue"],[37380,"ExposureCompensation"],[37381,"MaxApertureValue"],[37382,"SubjectDistance"],[37383,"MeteringMode"],[37384,"LightSource"],[37385,"Flash"],[37386,"FocalLength"],[37393,"ImageNumber"],[37394,"SecurityClassification"],[37395,"ImageHistory"],[37396,"SubjectArea"],[37500,"MakerNote"],[37510,"UserComment"],[37520,"SubSecTime"],[37521,"SubSecTimeOriginal"],[37522,"SubSecTimeDigitized"],[37888,"AmbientTemperature"],[37889,"Humidity"],[37890,"Pressure"],[37891,"WaterDepth"],[37892,"Acceleration"],[37893,"CameraElevationAngle"],[40960,"FlashpixVersion"],[40961,"ColorSpace"],[40962,"ExifImageWidth"],[40963,"ExifImageHeight"],[40964,"RelatedSoundFile"],[41483,"FlashEnergy"],[41486,"FocalPlaneXResolution"],[41487,"FocalPlaneYResolution"],[41488,"FocalPlaneResolutionUnit"],[41492,"SubjectLocation"],[41493,"ExposureIndex"],[41495,"SensingMethod"],[41728,"FileSource"],[41729,"SceneType"],[41730,"CFAPattern"],[41985,"CustomRendered"],[41986,"ExposureMode"],[41987,"WhiteBalance"],[41988,"DigitalZoomRatio"],[41989,"FocalLengthIn35mmFormat"],[41990,"SceneCaptureType"],[41991,"GainControl"],[41992,"Contrast"],[41993,"Saturation"],[41994,"Sharpness"],[41996,"SubjectDistanceRange"],[42016,"ImageUniqueID"],[42032,"OwnerName"],[42033,"SerialNumber"],[42034,"LensInfo"],[42035,"LensMake"],[42036,"LensModel"],[42037,"LensSerialNumber"],[42080,"CompositeImage"],[42081,"CompositeImageCount"],[42082,"CompositeImageExposureTimes"],[42240,"Gamma"],[59932,"Padding"],[59933,"OffsetSchema"],[65e3,"OwnerName"],[65001,"SerialNumber"],[65002,"Lens"],[65100,"RawFile"],[65101,"Converter"],[65102,"WhiteBalance"],[65105,"Exposure"],[65106,"Shadows"],[65107,"Brightness"],[65108,"Contrast"],[65109,"Saturation"],[65110,"Sharpness"],[65111,"Smoothness"],[65112,"MoireFilter"],[40965,"InteropIFD"]]),L(N,"gps",[[0,"GPSVersionID"],[1,"GPSLatitudeRef"],[2,"GPSLatitude"],[3,"GPSLongitudeRef"],[4,"GPSLongitude"],[5,"GPSAltitudeRef"],[6,"GPSAltitude"],[7,"GPSTimeStamp"],[8,"GPSSatellites"],[9,"GPSStatus"],[10,"GPSMeasureMode"],[11,"GPSDOP"],[12,"GPSSpeedRef"],[13,"GPSSpeed"],[14,"GPSTrackRef"],[15,"GPSTrack"],[16,"GPSImgDirectionRef"],[17,"GPSImgDirection"],[18,"GPSMapDatum"],[19,"GPSDestLatitudeRef"],[20,"GPSDestLatitude"],[21,"GPSDestLongitudeRef"],[22,"GPSDestLongitude"],[23,"GPSDestBearingRef"],[24,"GPSDestBearing"],[25,"GPSDestDistanceRef"],[26,"GPSDestDistance"],[27,"GPSProcessingMethod"],[28,"GPSAreaInformation"],[29,"GPSDateStamp"],[30,"GPSDifferential"],[31,"GPSHPositioningError"]]),L($,["ifd0","ifd1"],[[274,{1:"Horizontal (normal)",2:"Mirror horizontal",3:"Rotate 180",4:"Mirror vertical",5:"Mirror horizontal and rotate 270 CW",6:"Rotate 90 CW",7:"Mirror horizontal and rotate 90 CW",8:"Rotate 270 CW"}],[296,{1:"None",2:"inches",3:"cm"}]]);let De=L($,"exif",[[34850,{0:"Not defined",1:"Manual",2:"Normal program",3:"Aperture priority",4:"Shutter priority",5:"Creative program",6:"Action program",7:"Portrait mode",8:"Landscape mode"}],[37121,{0:"-",1:"Y",2:"Cb",3:"Cr",4:"R",5:"G",6:"B"}],[37383,{0:"Unknown",1:"Average",2:"CenterWeightedAverage",3:"Spot",4:"MultiSpot",5:"Pattern",6:"Partial",255:"Other"}],[37384,{0:"Unknown",1:"Daylight",2:"Fluorescent",3:"Tungsten (incandescent light)",4:"Flash",9:"Fine weather",10:"Cloudy weather",11:"Shade",12:"Daylight fluorescent (D 5700 - 7100K)",13:"Day white fluorescent (N 4600 - 5400K)",14:"Cool white fluorescent (W 3900 - 4500K)",15:"White fluorescent (WW 3200 - 3700K)",17:"Standard light A",18:"Standard light B",19:"Standard light C",20:"D55",21:"D65",22:"D75",23:"D50",24:"ISO studio tungsten",255:"Other"}],[37385,{0:"Flash did not fire",1:"Flash fired",5:"Strobe return light not detected",7:"Strobe return light detected",9:"Flash fired, compulsory flash mode",13:"Flash fired, compulsory flash mode, return light not detected",15:"Flash fired, compulsory flash mode, return light detected",16:"Flash did not fire, compulsory flash mode",24:"Flash did not fire, auto mode",25:"Flash fired, auto mode",29:"Flash fired, auto mode, return light not detected",31:"Flash fired, auto mode, return light detected",32:"No flash function",65:"Flash fired, red-eye reduction mode",69:"Flash fired, red-eye reduction mode, return light not detected",71:"Flash fired, red-eye reduction mode, return light detected",73:"Flash fired, compulsory flash mode, red-eye reduction mode",77:"Flash fired, compulsory flash mode, red-eye reduction mode, return light not detected",79:"Flash fired, compulsory flash mode, red-eye reduction mode, return light detected",89:"Flash fired, auto mode, red-eye reduction mode",93:"Flash fired, auto mode, return light not detected, red-eye reduction mode",95:"Flash fired, auto mode, return light detected, red-eye reduction mode"}],[41495,{1:"Not defined",2:"One-chip color area sensor",3:"Two-chip color area sensor",4:"Three-chip color area sensor",5:"Color sequential area sensor",7:"Trilinear sensor",8:"Color sequential linear sensor"}],[41728,{1:"Film Scanner",2:"Reflection Print Scanner",3:"Digital Camera"}],[41729,{1:"Directly photographed"}],[41985,{0:"Normal",1:"Custom",2:"HDR (no original saved)",3:"HDR (original saved)",4:"Original (for HDR)",6:"Panorama",7:"Portrait HDR",8:"Portrait"}],[41986,{0:"Auto",1:"Manual",2:"Auto bracket"}],[41987,{0:"Auto",1:"Manual"}],[41990,{0:"Standard",1:"Landscape",2:"Portrait",3:"Night",4:"Other"}],[41991,{0:"None",1:"Low gain up",2:"High gain up",3:"Low gain down",4:"High gain down"}],[41996,{0:"Unknown",1:"Macro",2:"Close",3:"Distant"}],[42080,{0:"Unknown",1:"Not a Composite Image",2:"General Composite Image",3:"Composite Image Captured While Shooting"}]]);const Ht={1:"No absolute unit of measurement",2:"Inch",3:"Centimeter"};De.set(37392,Ht),De.set(41488,Ht);const ht={0:"Normal",1:"Low",2:"High"};function jt(n){return typeof n=="object"&&n.length!==void 0?n[0]:n}function Wt(n){let e=Array.from(n).slice(1);return e[1]>15&&(e=e.map(t=>String.fromCharCode(t))),e[2]!=="0"&&e[2]!==0||e.pop(),e.join(".")}function pt(n){if(typeof n=="string"){var[e,t,r,i,s,a]=n.trim().split(/[-: ]/g).map(Number),o=new Date(e,t-1,r);return Number.isNaN(i)||Number.isNaN(s)||Number.isNaN(a)||(o.setHours(i),o.setMinutes(s),o.setSeconds(a)),Number.isNaN(+o)?n:o}}function Oe(n){if(typeof n=="string")return n;let e=[];if(n[1]===0&&n[n.length-1]===0)for(let t=0;t<n.length;t+=2)e.push(Xt(n[t+1],n[t]));else for(let t=0;t<n.length;t+=2)e.push(Xt(n[t],n[t+1]));return _e(String.fromCodePoint(...e))}function Xt(n,e){return n<<8|e}De.set(41992,ht),De.set(41993,ht),De.set(41994,ht),L(he,["ifd0","ifd1"],[[50827,function(n){return typeof n!="string"?ar(n):n}],[306,pt],[40091,Oe],[40092,Oe],[40093,Oe],[40094,Oe],[40095,Oe]]),L(he,"exif",[[40960,Wt],[36864,Wt],[36867,pt],[36868,pt],[40962,jt],[40963,jt]]),L(he,"gps",[[0,n=>Array.from(n).join(".")],[7,n=>Array.from(n).join(":")]]);class ut extends Y{static canHandle(e,t){return e.getUint8(t+1)===225&&e.getUint32(t+4)===1752462448&&e.getString(t+4,20)==="http://ns.adobe.com/"}static headerLength(e,t){return e.getString(t+4,34)==="http://ns.adobe.com/xmp/extension/"?79:33}static findPosition(e,t){let r=super.findPosition(e,t);return r.multiSegment=r.extended=r.headerLength===79,r.multiSegment?(r.chunkCount=e.getUint8(t+72),r.chunkNumber=e.getUint8(t+76),e.getUint8(t+77)!==0&&r.chunkNumber++):(r.chunkCount=1/0,r.chunkNumber=-1),r}static handleMultiSegments(e){return e.map(t=>t.chunk.getString()).join("")}normalizeInput(e){return typeof e=="string"?e:j.from(e).getString()}parse(e=this.chunk){if(!this.localOptions.parse)return e;e=function(s){let a={},o={};for(let l of _r)a[l]=[],o[l]=0;return s.replace(dn,(l,h,c)=>{if(h==="<"){let d=++o[c];return a[c].push(d),`${l}#${d}`}return`${l}#${a[c].pop()}`})}(e);let t=Pe.findAll(e,"rdf","Description");t.length===0&&t.push(new Pe("rdf","Description",void 0,e));let r,i={};for(let s of t)for(let a of s.properties)r=ln(a.ns,i),fr(a,r);return function(s){let a;for(let o in s)a=s[o]=Xe(s[o]),a===void 0&&delete s[o];return Xe(s)}(i)}assignToOutput(e,t){if(this.localOptions.parse)for(let[r,i]of Object.entries(t))switch(r){case"tiff":this.assignObjectToOutput(e,"ifd0",i);break;case"exif":this.assignObjectToOutput(e,"exif",i);break;case"xmlns":break;default:this.assignObjectToOutput(e,r,i)}else e.xmp=t}}m(ut,"type","xmp"),m(ut,"multiSegment",!0),F.set("xmp",ut);class qe{static findAll(e){return gr(e,/([a-zA-Z0-9-]+):([a-zA-Z0-9-]+)=("[^"]*"|'[^']*')/gm).map(qe.unpackMatch)}static unpackMatch(e){let t=e[1],r=e[2],i=e[3].slice(1,-1);return i=mr(i),new qe(t,r,i)}constructor(e,t,r){this.ns=e,this.name=t,this.value=r}serialize(){return this.value}}class Pe{static findAll(e,t,r){if(t!==void 0||r!==void 0){t=t||"[\\w\\d-]+",r=r||"[\\w\\d-]+";var i=new RegExp(`<(${t}):(${r})(#\\d+)?((\\s+?[\\w\\d-:]+=("[^"]*"|'[^']*'))*\\s*)(\\/>|>([\\s\\S]*?)<\\/\\1:\\2\\3>)`,"gm")}else i=/<([\w\d-]+):([\w\d-]+)(#\d+)?((\s+?[\w\d-:]+=("[^"]*"|'[^']*'))*\s*)(\/>|>([\s\S]*?)<\/\1:\2\3>)/gm;return gr(e,i).map(Pe.unpackMatch)}static unpackMatch(e){let t=e[1],r=e[2],i=e[4],s=e[8];return new Pe(t,r,i,s)}constructor(e,t,r,i){this.ns=e,this.name=t,this.attrString=r,this.innerXml=i,this.attrs=qe.findAll(r),this.children=Pe.findAll(i),this.value=this.children.length===0?mr(i):void 0,this.properties=[...this.attrs,...this.children]}get isPrimitive(){return this.value!==void 0&&this.attrs.length===0&&this.children.length===0}get isListContainer(){return this.children.length===1&&this.children[0].isList}get isList(){let{ns:e,name:t}=this;return e==="rdf"&&(t==="Seq"||t==="Bag"||t==="Alt")}get isListItem(){return this.ns==="rdf"&&this.name==="li"}serialize(){if(this.properties.length===0&&this.value===void 0)return;if(this.isPrimitive)return this.value;if(this.isListContainer)return this.children[0].serialize();if(this.isList)return on(this.children.map(an));if(this.isListItem&&this.children.length===1&&this.attrs.length===0)return this.children[0].serialize();let e={};for(let t of this.properties)fr(t,e);return this.value!==void 0&&(e.value=this.value),Xe(e)}}function fr(n,e){let t=n.serialize();t!==void 0&&(e[n.name]=t)}var an=n=>n.serialize(),on=n=>n.length===1?n[0]:n,ln=(n,e)=>e[n]?e[n]:e[n]={};function gr(n,e){let t,r=[];if(!n)return r;for(;(t=e.exec(n))!==null;)r.push(t);return r}function mr(n){if(function(r){return r==null||r==="null"||r==="undefined"||r===""||r.trim()===""}(n))return;let e=Number(n);if(!Number.isNaN(e))return e;let t=n.toLowerCase();return t==="true"||t!=="false"&&n.trim()}const _r=["rdf:li","rdf:Seq","rdf:Bag","rdf:Alt","rdf:Description"],dn=new RegExp(`(<|\\/)(${_r.join("|")})`,"g");var yr=Object.freeze({__proto__:null,default:nn,Exifr:pe,fileParsers:ee,segmentParsers:F,fileReaders:te,tagKeys:N,tagValues:$,tagRevivers:he,createDictionary:L,extendDictionary:Se,fetchUrlAsArrayBuffer:be,readBlobAsArrayBuffer:xe,chunkedProps:oe,otherSegments:Ae,segments:ke,tiffBlocks:I,segmentsAndBlocks:le,tiffExtractables:de,inheritables:Ie,allFormatters:ce,Options:Ce,parse:Ze,gpsOnlyOptions:Ct,gps:lr,thumbnailOnlyOptions:vt,thumbnail:dr,thumbnailUrl:cr,orientationOnlyOptions:Pt,orientation:At,rotations:It,get rotateCanvas(){return ye},get rotateCss(){return we},rotation:hr});let Kt=We("fs",n=>n.promises);te.set("fs",class extends tt{async readWhole(){this.chunked=!1,this.fs=await Kt;let n=await this.fs.readFile(this.input);this._swapBuffer(n)}async readChunked(){this.chunked=!0,this.fs=await Kt,await this.open(),await this.readChunk(0,this.options.firstChunkSize)}async open(){this.fh===void 0&&(this.fh=await this.fs.open(this.input,"r"),this.size=(await this.fh.stat(this.input)).size)}async _readChunk(n,e){this.fh===void 0&&await this.open(),n+e>this.size&&(e=this.size-n);var t=this.subarray(n,e,!0);return await this.fh.read(t.dataView,0,e,n),t}async close(){if(this.fh){let n=this.fh;this.fh=void 0,await n.close()}}});te.set("base64",class extends tt{constructor(...n){super(...n),this.input=this.input.replace(/^data:([^;]+);base64,/gim,""),this.size=this.input.length/4*3,this.input.endsWith("==")?this.size-=2:this.input.endsWith("=")&&(this.size-=1)}async _readChunk(n,e){let t,r,i=this.input;n===void 0?(n=0,t=0,r=0):(t=4*Math.floor(n/3),r=n-t/4*3),e===void 0&&(e=this.size);let s=n+e,a=t+4*Math.ceil(s/3);i=i.slice(t,a);let o=Math.min(e,this.size-n);if(Qe){let l=Je.from(i,"base64").slice(r,r+o);return this.set(l,n,!0)}{let l=this.subarray(n,o,!0),h=atob(i),c=l.toUint8();for(let d=0;d<o;d++)c[d]=h.charCodeAt(r+d);return l}}});class qt extends et{static canHandle(e,t){return t===18761||t===19789}extendOptions(e){let{ifd0:t,xmp:r,iptc:i,icc:s}=e;r.enabled&&t.deps.add(700),i.enabled&&t.deps.add(33723),s.enabled&&t.deps.add(34675),t.finalizeFilters()}async parse(){let{tiff:e,xmp:t,iptc:r,icc:i}=this.options;if(e.enabled||t.enabled||r.enabled||i.enabled){let s=Math.max(gt(this.options),this.options.chunkSize);await this.file.ensureChunk(0,s),this.createParser("tiff",this.file),this.parsers.tiff.parseHeader(),await this.parsers.tiff.parseIfd0Block(),this.adaptTiffPropAsSegment("xmp"),this.adaptTiffPropAsSegment("iptc"),this.adaptTiffPropAsSegment("icc")}}adaptTiffPropAsSegment(e){if(this.parsers.tiff[e]){let t=this.parsers.tiff[e];this.injectSegment(e,t)}}}m(qt,"type","tiff"),ee.set("tiff",qt);let cn=We("zlib");const hn=["ihdr","iccp","text","itxt","exif"];class Yt extends et{constructor(...e){super(...e),m(this,"catchError",t=>this.errors.push(t)),m(this,"metaChunks",[]),m(this,"unknownChunks",[])}static canHandle(e,t){return t===35152&&e.getUint32(0)===2303741511&&e.getUint32(4)===218765834}async parse(){let{file:e}=this;await this.findPngChunksInRange(8,e.byteLength),await this.readSegments(this.metaChunks),this.findIhdr(),this.parseTextChunks(),await this.findExif().catch(this.catchError),await this.findXmp().catch(this.catchError),await this.findIcc().catch(this.catchError)}async findPngChunksInRange(e,t){let{file:r}=this;for(;e<t;){let i=r.getUint32(e),s=r.getUint32(e+4),a=r.getString(e+4,4).toLowerCase(),o=i+4+4+4,l={type:a,offset:e,length:o,start:e+4+4,size:i,marker:s};hn.includes(a)?this.metaChunks.push(l):this.unknownChunks.push(l),e+=o}}parseTextChunks(){let e=this.metaChunks.filter(t=>t.type==="text");for(let t of e){let[r,i]=this.file.getString(t.start,t.size).split("\0");this.injectKeyValToIhdr(r,i)}}injectKeyValToIhdr(e,t){let r=this.parsers.ihdr;r&&r.raw.set(e,t)}findIhdr(){let e=this.metaChunks.find(t=>t.type==="ihdr");e&&this.options.ihdr.enabled!==!1&&this.createParser("ihdr",e.chunk)}async findExif(){let e=this.metaChunks.find(t=>t.type==="exif");e&&this.injectSegment("tiff",e.chunk)}async findXmp(){let e=this.metaChunks.filter(t=>t.type==="itxt");for(let t of e)t.chunk.getString(0,17)==="XML:com.adobe.xmp"&&this.injectSegment("xmp",t.chunk)}async findIcc(){let e=this.metaChunks.find(o=>o.type==="iccp");if(!e)return;let{chunk:t}=e,r=t.getUint8Array(0,81),i=0;for(;i<80&&r[i]!==0;)i++;let s=i+2,a=t.getString(0,i);if(this.injectKeyValToIhdr("ProfileName",a),je){let o=await cn,l=t.getUint8Array(s);l=o.inflateSync(l),this.injectSegment("icc",l)}}}m(Yt,"type","png"),ee.set("png",Yt),L(N,"interop",[[1,"InteropIndex"],[2,"InteropVersion"],[4096,"RelatedImageFileFormat"],[4097,"RelatedImageWidth"],[4098,"RelatedImageHeight"]]),Se(N,"ifd0",[[11,"ProcessingSoftware"],[254,"SubfileType"],[255,"OldSubfileType"],[263,"Thresholding"],[264,"CellWidth"],[265,"CellLength"],[266,"FillOrder"],[269,"DocumentName"],[280,"MinSampleValue"],[281,"MaxSampleValue"],[285,"PageName"],[286,"XPosition"],[287,"YPosition"],[290,"GrayResponseUnit"],[297,"PageNumber"],[321,"HalftoneHints"],[322,"TileWidth"],[323,"TileLength"],[332,"InkSet"],[337,"TargetPrinter"],[18246,"Rating"],[18249,"RatingPercent"],[33550,"PixelScale"],[34264,"ModelTransform"],[34377,"PhotoshopSettings"],[50706,"DNGVersion"],[50707,"DNGBackwardVersion"],[50708,"UniqueCameraModel"],[50709,"LocalizedCameraModel"],[50736,"DNGLensInfo"],[50739,"ShadowScale"],[50740,"DNGPrivateData"],[33920,"IntergraphMatrix"],[33922,"ModelTiePoint"],[34118,"SEMInfo"],[34735,"GeoTiffDirectory"],[34736,"GeoTiffDoubleParams"],[34737,"GeoTiffAsciiParams"],[50341,"PrintIM"],[50721,"ColorMatrix1"],[50722,"ColorMatrix2"],[50723,"CameraCalibration1"],[50724,"CameraCalibration2"],[50725,"ReductionMatrix1"],[50726,"ReductionMatrix2"],[50727,"AnalogBalance"],[50728,"AsShotNeutral"],[50729,"AsShotWhiteXY"],[50730,"BaselineExposure"],[50731,"BaselineNoise"],[50732,"BaselineSharpness"],[50734,"LinearResponseLimit"],[50735,"CameraSerialNumber"],[50741,"MakerNoteSafety"],[50778,"CalibrationIlluminant1"],[50779,"CalibrationIlluminant2"],[50781,"RawDataUniqueID"],[50827,"OriginalRawFileName"],[50828,"OriginalRawFileData"],[50831,"AsShotICCProfile"],[50832,"AsShotPreProfileMatrix"],[50833,"CurrentICCProfile"],[50834,"CurrentPreProfileMatrix"],[50879,"ColorimetricReference"],[50885,"SRawType"],[50898,"PanasonicTitle"],[50899,"PanasonicTitle2"],[50931,"CameraCalibrationSig"],[50932,"ProfileCalibrationSig"],[50933,"ProfileIFD"],[50934,"AsShotProfileName"],[50936,"ProfileName"],[50937,"ProfileHueSatMapDims"],[50938,"ProfileHueSatMapData1"],[50939,"ProfileHueSatMapData2"],[50940,"ProfileToneCurve"],[50941,"ProfileEmbedPolicy"],[50942,"ProfileCopyright"],[50964,"ForwardMatrix1"],[50965,"ForwardMatrix2"],[50966,"PreviewApplicationName"],[50967,"PreviewApplicationVersion"],[50968,"PreviewSettingsName"],[50969,"PreviewSettingsDigest"],[50970,"PreviewColorSpace"],[50971,"PreviewDateTime"],[50972,"RawImageDigest"],[50973,"OriginalRawFileDigest"],[50981,"ProfileLookTableDims"],[50982,"ProfileLookTableData"],[51043,"TimeCodes"],[51044,"FrameRate"],[51058,"TStop"],[51081,"ReelName"],[51089,"OriginalDefaultFinalSize"],[51090,"OriginalBestQualitySize"],[51091,"OriginalDefaultCropSize"],[51105,"CameraLabel"],[51107,"ProfileHueSatMapEncoding"],[51108,"ProfileLookTableEncoding"],[51109,"BaselineExposureOffset"],[51110,"DefaultBlackRender"],[51111,"NewRawImageDigest"],[51112,"RawToPreviewGain"]]);let $t=[[273,"StripOffsets"],[279,"StripByteCounts"],[288,"FreeOffsets"],[289,"FreeByteCounts"],[291,"GrayResponseCurve"],[292,"T4Options"],[293,"T6Options"],[300,"ColorResponseUnit"],[320,"ColorMap"],[324,"TileOffsets"],[325,"TileByteCounts"],[326,"BadFaxLines"],[327,"CleanFaxData"],[328,"ConsecutiveBadFaxLines"],[330,"SubIFD"],[333,"InkNames"],[334,"NumberofInks"],[336,"DotRange"],[338,"ExtraSamples"],[339,"SampleFormat"],[340,"SMinSampleValue"],[341,"SMaxSampleValue"],[342,"TransferRange"],[343,"ClipPath"],[344,"XClipPathUnits"],[345,"YClipPathUnits"],[346,"Indexed"],[347,"JPEGTables"],[351,"OPIProxy"],[400,"GlobalParametersIFD"],[401,"ProfileType"],[402,"FaxProfile"],[403,"CodingMethods"],[404,"VersionYear"],[405,"ModeNumber"],[433,"Decode"],[434,"DefaultImageColor"],[435,"T82Options"],[437,"JPEGTables"],[512,"JPEGProc"],[515,"JPEGRestartInterval"],[517,"JPEGLosslessPredictors"],[518,"JPEGPointTransforms"],[519,"JPEGQTables"],[520,"JPEGDCTables"],[521,"JPEGACTables"],[559,"StripRowCounts"],[999,"USPTOMiscellaneous"],[18247,"XP_DIP_XML"],[18248,"StitchInfo"],[28672,"SonyRawFileType"],[28688,"SonyToneCurve"],[28721,"VignettingCorrection"],[28722,"VignettingCorrParams"],[28724,"ChromaticAberrationCorrection"],[28725,"ChromaticAberrationCorrParams"],[28726,"DistortionCorrection"],[28727,"DistortionCorrParams"],[29895,"SonyCropTopLeft"],[29896,"SonyCropSize"],[32781,"ImageID"],[32931,"WangTag1"],[32932,"WangAnnotation"],[32933,"WangTag3"],[32934,"WangTag4"],[32953,"ImageReferencePoints"],[32954,"RegionXformTackPoint"],[32955,"WarpQuadrilateral"],[32956,"AffineTransformMat"],[32995,"Matteing"],[32996,"DataType"],[32997,"ImageDepth"],[32998,"TileDepth"],[33300,"ImageFullWidth"],[33301,"ImageFullHeight"],[33302,"TextureFormat"],[33303,"WrapModes"],[33304,"FovCot"],[33305,"MatrixWorldToScreen"],[33306,"MatrixWorldToCamera"],[33405,"Model2"],[33421,"CFARepeatPatternDim"],[33422,"CFAPattern2"],[33423,"BatteryLevel"],[33424,"KodakIFD"],[33445,"MDFileTag"],[33446,"MDScalePixel"],[33447,"MDColorTable"],[33448,"MDLabName"],[33449,"MDSampleInfo"],[33450,"MDPrepDate"],[33451,"MDPrepTime"],[33452,"MDFileUnits"],[33589,"AdventScale"],[33590,"AdventRevision"],[33628,"UIC1Tag"],[33629,"UIC2Tag"],[33630,"UIC3Tag"],[33631,"UIC4Tag"],[33918,"IntergraphPacketData"],[33919,"IntergraphFlagRegisters"],[33921,"INGRReserved"],[34016,"Site"],[34017,"ColorSequence"],[34018,"IT8Header"],[34019,"RasterPadding"],[34020,"BitsPerRunLength"],[34021,"BitsPerExtendedRunLength"],[34022,"ColorTable"],[34023,"ImageColorIndicator"],[34024,"BackgroundColorIndicator"],[34025,"ImageColorValue"],[34026,"BackgroundColorValue"],[34027,"PixelIntensityRange"],[34028,"TransparencyIndicator"],[34029,"ColorCharacterization"],[34030,"HCUsage"],[34031,"TrapIndicator"],[34032,"CMYKEquivalent"],[34152,"AFCP_IPTC"],[34232,"PixelMagicJBIGOptions"],[34263,"JPLCartoIFD"],[34306,"WB_GRGBLevels"],[34310,"LeafData"],[34687,"TIFF_FXExtensions"],[34688,"MultiProfiles"],[34689,"SharedData"],[34690,"T88Options"],[34732,"ImageLayer"],[34750,"JBIGOptions"],[34856,"Opto-ElectricConvFactor"],[34857,"Interlace"],[34908,"FaxRecvParams"],[34909,"FaxSubAddress"],[34910,"FaxRecvTime"],[34929,"FedexEDR"],[34954,"LeafSubIFD"],[37387,"FlashEnergy"],[37388,"SpatialFrequencyResponse"],[37389,"Noise"],[37390,"FocalPlaneXResolution"],[37391,"FocalPlaneYResolution"],[37392,"FocalPlaneResolutionUnit"],[37397,"ExposureIndex"],[37398,"TIFF-EPStandardID"],[37399,"SensingMethod"],[37434,"CIP3DataFile"],[37435,"CIP3Sheet"],[37436,"CIP3Side"],[37439,"StoNits"],[37679,"MSDocumentText"],[37680,"MSPropertySetStorage"],[37681,"MSDocumentTextPosition"],[37724,"ImageSourceData"],[40965,"InteropIFD"],[40976,"SamsungRawPointersOffset"],[40977,"SamsungRawPointersLength"],[41217,"SamsungRawByteOrder"],[41218,"SamsungRawUnknown"],[41484,"SpatialFrequencyResponse"],[41485,"Noise"],[41489,"ImageNumber"],[41490,"SecurityClassification"],[41491,"ImageHistory"],[41494,"TIFF-EPStandardID"],[41995,"DeviceSettingDescription"],[42112,"GDALMetadata"],[42113,"GDALNoData"],[44992,"ExpandSoftware"],[44993,"ExpandLens"],[44994,"ExpandFilm"],[44995,"ExpandFilterLens"],[44996,"ExpandScanner"],[44997,"ExpandFlashLamp"],[46275,"HasselbladRawImage"],[48129,"PixelFormat"],[48130,"Transformation"],[48131,"Uncompressed"],[48132,"ImageType"],[48256,"ImageWidth"],[48257,"ImageHeight"],[48258,"WidthResolution"],[48259,"HeightResolution"],[48320,"ImageOffset"],[48321,"ImageByteCount"],[48322,"AlphaOffset"],[48323,"AlphaByteCount"],[48324,"ImageDataDiscard"],[48325,"AlphaDataDiscard"],[50215,"OceScanjobDesc"],[50216,"OceApplicationSelector"],[50217,"OceIDNumber"],[50218,"OceImageLogic"],[50255,"Annotations"],[50459,"HasselbladExif"],[50547,"OriginalFileName"],[50560,"USPTOOriginalContentType"],[50656,"CR2CFAPattern"],[50710,"CFAPlaneColor"],[50711,"CFALayout"],[50712,"LinearizationTable"],[50713,"BlackLevelRepeatDim"],[50714,"BlackLevel"],[50715,"BlackLevelDeltaH"],[50716,"BlackLevelDeltaV"],[50717,"WhiteLevel"],[50718,"DefaultScale"],[50719,"DefaultCropOrigin"],[50720,"DefaultCropSize"],[50733,"BayerGreenSplit"],[50737,"ChromaBlurRadius"],[50738,"AntiAliasStrength"],[50752,"RawImageSegmentation"],[50780,"BestQualityScale"],[50784,"AliasLayerMetadata"],[50829,"ActiveArea"],[50830,"MaskedAreas"],[50935,"NoiseReductionApplied"],[50974,"SubTileBlockSize"],[50975,"RowInterleaveFactor"],[51008,"OpcodeList1"],[51009,"OpcodeList2"],[51022,"OpcodeList3"],[51041,"NoiseProfile"],[51114,"CacheVersion"],[51125,"DefaultUserCrop"],[51157,"NikonNEFInfo"],[65024,"KdcIFD"]];Se(N,"ifd0",$t),Se(N,"exif",$t),L($,"gps",[[23,{M:"Magnetic North",T:"True North"}],[25,{K:"Kilometers",M:"Miles",N:"Nautical Miles"}]]);class ft extends Y{static canHandle(e,t){return e.getUint8(t+1)===224&&e.getUint32(t+4)===1246120262&&e.getUint8(t+8)===0}parse(){return this.parseTags(),this.translate(),this.output}parseTags(){this.raw=new Map([[0,this.chunk.getUint16(0)],[2,this.chunk.getUint8(2)],[3,this.chunk.getUint16(3)],[5,this.chunk.getUint16(5)],[7,this.chunk.getUint8(7)],[8,this.chunk.getUint8(8)]])}}m(ft,"type","jfif"),m(ft,"headerLength",9),F.set("jfif",ft),L(N,"jfif",[[0,"JFIFVersion"],[2,"ResolutionUnit"],[3,"XResolution"],[5,"YResolution"],[7,"ThumbnailWidth"],[8,"ThumbnailHeight"]]);class Jt extends Y{parse(){return this.parseTags(),this.translate(),this.output}parseTags(){this.raw=new Map([[0,this.chunk.getUint32(0)],[4,this.chunk.getUint32(4)],[8,this.chunk.getUint8(8)],[9,this.chunk.getUint8(9)],[10,this.chunk.getUint8(10)],[11,this.chunk.getUint8(11)],[12,this.chunk.getUint8(12)],...Array.from(this.raw)])}}m(Jt,"type","ihdr"),F.set("ihdr",Jt),L(N,"ihdr",[[0,"ImageWidth"],[4,"ImageHeight"],[8,"BitDepth"],[9,"ColorType"],[10,"Compression"],[11,"Filter"],[12,"Interlace"]]),L($,"ihdr",[[9,{0:"Grayscale",2:"RGB",3:"Palette",4:"Grayscale with Alpha",6:"RGB with Alpha",DEFAULT:"Unknown"}],[10,{0:"Deflate/Inflate",DEFAULT:"Unknown"}],[11,{0:"Adaptive",DEFAULT:"Unknown"}],[12,{0:"Noninterlaced",1:"Adam7 Interlace",DEFAULT:"Unknown"}]]);class ze extends Y{static canHandle(e,t){return e.getUint8(t+1)===226&&e.getUint32(t+4)===1229144927}static findPosition(e,t){let r=super.findPosition(e,t);return r.chunkNumber=e.getUint8(t+16),r.chunkCount=e.getUint8(t+17),r.multiSegment=r.chunkCount>1,r}static handleMultiSegments(e){return function(t){let r=function(i){let s=i[0].constructor,a=0;for(let h of i)a+=h.length;let o=new s(a),l=0;for(let h of i)o.set(h,l),l+=h.length;return o}(t.map(i=>i.chunk.toUint8()));return new j(r)}(e)}parse(){return this.raw=new Map,this.parseHeader(),this.parseTags(),this.translate(),this.output}parseHeader(){let{raw:e}=this;this.chunk.byteLength<84&&M("ICC header is too short");for(let[t,r]of Object.entries(pn)){t=parseInt(t,10);let i=r(this.chunk,t);i!=="\0\0\0\0"&&e.set(t,i)}}parseTags(){let e,t,r,i,s,{raw:a}=this,o=this.chunk.getUint32(128),l=132,h=this.chunk.byteLength;for(;o--;){if(e=this.chunk.getString(l,4),t=this.chunk.getUint32(l+4),r=this.chunk.getUint32(l+8),i=this.chunk.getString(t,4),t+r>h)return void console.warn("reached the end of the first ICC chunk. Enable options.tiff.multiSegment to read all ICC segments.");s=this.parseTag(i,t,r),s!==void 0&&s!=="\0\0\0\0"&&a.set(e,s),l+=12}}parseTag(e,t,r){switch(e){case"desc":return this.parseDesc(t);case"mluc":return this.parseMluc(t);case"text":return this.parseText(t,r);case"sig ":return this.parseSig(t)}if(!(t+r>this.chunk.byteLength))return this.chunk.getUint8Array(t,r)}parseDesc(e){let t=this.chunk.getUint32(e+8)-1;return _e(this.chunk.getString(e+12,t))}parseText(e,t){return _e(this.chunk.getString(e+8,t-8))}parseSig(e){return _e(this.chunk.getString(e+8,4))}parseMluc(e){let{chunk:t}=this,r=t.getUint32(e+8),i=t.getUint32(e+12),s=e+16,a=[];for(let o=0;o<r;o++){let l=t.getString(s+0,2),h=t.getString(s+2,2),c=t.getUint32(s+4),d=t.getUint32(s+8)+e,u=_e(t.getUnicodeString(d,c));a.push({lang:l,country:h,text:u}),s+=i}return r===1?a[0].text:a}translateValue(e,t){return typeof e=="string"?t[e]||t[e.toLowerCase()]||e:t[e]||e}}m(ze,"type","icc"),m(ze,"multiSegment",!0),m(ze,"headerLength",18);const pn={4:ie,8:function(n,e){return[n.getUint8(e),n.getUint8(e+1)>>4,n.getUint8(e+1)%16].map(t=>t.toString(10)).join(".")},12:ie,16:ie,20:ie,24:function(n,e){const t=n.getUint16(e),r=n.getUint16(e+2)-1,i=n.getUint16(e+4),s=n.getUint16(e+6),a=n.getUint16(e+8),o=n.getUint16(e+10);return new Date(Date.UTC(t,r,i,s,a,o))},36:ie,40:ie,48:ie,52:ie,64:(n,e)=>n.getUint32(e),80:ie};function ie(n,e){return _e(n.getString(e,4))}F.set("icc",ze),L(N,"icc",[[4,"ProfileCMMType"],[8,"ProfileVersion"],[12,"ProfileClass"],[16,"ColorSpaceData"],[20,"ProfileConnectionSpace"],[24,"ProfileDateTime"],[36,"ProfileFileSignature"],[40,"PrimaryPlatform"],[44,"CMMFlags"],[48,"DeviceManufacturer"],[52,"DeviceModel"],[56,"DeviceAttributes"],[64,"RenderingIntent"],[68,"ConnectionSpaceIlluminant"],[80,"ProfileCreator"],[84,"ProfileID"],["Header","ProfileHeader"],["MS00","WCSProfiles"],["bTRC","BlueTRC"],["bXYZ","BlueMatrixColumn"],["bfd","UCRBG"],["bkpt","MediaBlackPoint"],["calt","CalibrationDateTime"],["chad","ChromaticAdaptation"],["chrm","Chromaticity"],["ciis","ColorimetricIntentImageState"],["clot","ColorantTableOut"],["clro","ColorantOrder"],["clrt","ColorantTable"],["cprt","ProfileCopyright"],["crdi","CRDInfo"],["desc","ProfileDescription"],["devs","DeviceSettings"],["dmdd","DeviceModelDesc"],["dmnd","DeviceMfgDesc"],["dscm","ProfileDescriptionML"],["fpce","FocalPlaneColorimetryEstimates"],["gTRC","GreenTRC"],["gXYZ","GreenMatrixColumn"],["gamt","Gamut"],["kTRC","GrayTRC"],["lumi","Luminance"],["meas","Measurement"],["meta","Metadata"],["mmod","MakeAndModel"],["ncl2","NamedColor2"],["ncol","NamedColor"],["ndin","NativeDisplayInfo"],["pre0","Preview0"],["pre1","Preview1"],["pre2","Preview2"],["ps2i","PS2RenderingIntent"],["ps2s","PostScript2CSA"],["psd0","PostScript2CRD0"],["psd1","PostScript2CRD1"],["psd2","PostScript2CRD2"],["psd3","PostScript2CRD3"],["pseq","ProfileSequenceDesc"],["psid","ProfileSequenceIdentifier"],["psvm","PS2CRDVMSize"],["rTRC","RedTRC"],["rXYZ","RedMatrixColumn"],["resp","OutputResponse"],["rhoc","ReflectionHardcopyOrigColorimetry"],["rig0","PerceptualRenderingIntentGamut"],["rig2","SaturationRenderingIntentGamut"],["rpoc","ReflectionPrintOutputColorimetry"],["sape","SceneAppearanceEstimates"],["scoe","SceneColorimetryEstimates"],["scrd","ScreeningDesc"],["scrn","Screening"],["targ","CharTarget"],["tech","Technology"],["vcgt","VideoCardGamma"],["view","ViewingConditions"],["vued","ViewingCondDesc"],["wtpt","MediaWhitePoint"]]);const Be={"4d2p":"Erdt Systems",AAMA:"Aamazing Technologies",ACER:"Acer",ACLT:"Acolyte Color Research",ACTI:"Actix Sytems",ADAR:"Adara Technology",ADBE:"Adobe",ADI:"ADI Systems",AGFA:"Agfa Graphics",ALMD:"Alps Electric",ALPS:"Alps Electric",ALWN:"Alwan Color Expertise",AMTI:"Amiable Technologies",AOC:"AOC International",APAG:"Apago",APPL:"Apple Computer",AST:"AST","AT&T":"AT&T",BAEL:"BARBIERI electronic",BRCO:"Barco NV",BRKP:"Breakpoint",BROT:"Brother",BULL:"Bull",BUS:"Bus Computer Systems","C-IT":"C-Itoh",CAMR:"Intel",CANO:"Canon",CARR:"Carroll Touch",CASI:"Casio",CBUS:"Colorbus PL",CEL:"Crossfield",CELx:"Crossfield",CGS:"CGS Publishing Technologies International",CHM:"Rochester Robotics",CIGL:"Colour Imaging Group, London",CITI:"Citizen",CL00:"Candela",CLIQ:"Color IQ",CMCO:"Chromaco",CMiX:"CHROMiX",COLO:"Colorgraphic Communications",COMP:"Compaq",COMp:"Compeq/Focus Technology",CONR:"Conrac Display Products",CORD:"Cordata Technologies",CPQ:"Compaq",CPRO:"ColorPro",CRN:"Cornerstone",CTX:"CTX International",CVIS:"ColorVision",CWC:"Fujitsu Laboratories",DARI:"Darius Technology",DATA:"Dataproducts",DCP:"Dry Creek Photo",DCRC:"Digital Contents Resource Center, Chung-Ang University",DELL:"Dell Computer",DIC:"Dainippon Ink and Chemicals",DICO:"Diconix",DIGI:"Digital","DL&C":"Digital Light & Color",DPLG:"Doppelganger",DS:"Dainippon Screen",DSOL:"DOOSOL",DUPN:"DuPont",EPSO:"Epson",ESKO:"Esko-Graphics",ETRI:"Electronics and Telecommunications Research Institute",EVER:"Everex Systems",EXAC:"ExactCODE",Eizo:"Eizo",FALC:"Falco Data Products",FF:"Fuji Photo Film",FFEI:"FujiFilm Electronic Imaging",FNRD:"Fnord Software",FORA:"Fora",FORE:"Forefront Technology",FP:"Fujitsu",FPA:"WayTech Development",FUJI:"Fujitsu",FX:"Fuji Xerox",GCC:"GCC Technologies",GGSL:"Global Graphics Software",GMB:"Gretagmacbeth",GMG:"GMG",GOLD:"GoldStar Technology",GOOG:"Google",GPRT:"Giantprint",GTMB:"Gretagmacbeth",GVC:"WayTech Development",GW2K:"Sony",HCI:"HCI",HDM:"Heidelberger Druckmaschinen",HERM:"Hermes",HITA:"Hitachi America",HP:"Hewlett-Packard",HTC:"Hitachi",HiTi:"HiTi Digital",IBM:"IBM",IDNT:"Scitex",IEC:"Hewlett-Packard",IIYA:"Iiyama North America",IKEG:"Ikegami Electronics",IMAG:"Image Systems",IMI:"Ingram Micro",INTC:"Intel",INTL:"N/A (INTL)",INTR:"Intra Electronics",IOCO:"Iocomm International Technology",IPS:"InfoPrint Solutions Company",IRIS:"Scitex",ISL:"Ichikawa Soft Laboratory",ITNL:"N/A (ITNL)",IVM:"IVM",IWAT:"Iwatsu Electric",Idnt:"Scitex",Inca:"Inca Digital Printers",Iris:"Scitex",JPEG:"Joint Photographic Experts Group",JSFT:"Jetsoft Development",JVC:"JVC Information Products",KART:"Scitex",KFC:"KFC Computek Components",KLH:"KLH Computers",KMHD:"Konica Minolta",KNCA:"Konica",KODA:"Kodak",KYOC:"Kyocera",Kart:"Scitex",LCAG:"Leica",LCCD:"Leeds Colour",LDAK:"Left Dakota",LEAD:"Leading Technology",LEXM:"Lexmark International",LINK:"Link Computer",LINO:"Linotronic",LITE:"Lite-On",Leaf:"Leaf",Lino:"Linotronic",MAGC:"Mag Computronic",MAGI:"MAG Innovision",MANN:"Mannesmann",MICN:"Micron Technology",MICR:"Microtek",MICV:"Microvitec",MINO:"Minolta",MITS:"Mitsubishi Electronics America",MITs:"Mitsuba",MNLT:"Minolta",MODG:"Modgraph",MONI:"Monitronix",MONS:"Monaco Systems",MORS:"Morse Technology",MOTI:"Motive Systems",MSFT:"Microsoft",MUTO:"MUTOH INDUSTRIES",Mits:"Mitsubishi Electric",NANA:"NANAO",NEC:"NEC",NEXP:"NexPress Solutions",NISS:"Nissei Sangyo America",NKON:"Nikon",NONE:"none",OCE:"Oce Technologies",OCEC:"OceColor",OKI:"Oki",OKID:"Okidata",OKIP:"Okidata",OLIV:"Olivetti",OLYM:"Olympus",ONYX:"Onyx Graphics",OPTI:"Optiquest",PACK:"Packard Bell",PANA:"Matsushita Electric Industrial",PANT:"Pantone",PBN:"Packard Bell",PFU:"PFU",PHIL:"Philips Consumer Electronics",PNTX:"HOYA",POne:"Phase One A/S",PREM:"Premier Computer Innovations",PRIN:"Princeton Graphic Systems",PRIP:"Princeton Publishing Labs",QLUX:"Hong Kong",QMS:"QMS",QPCD:"QPcard AB",QUAD:"QuadLaser",QUME:"Qume",RADI:"Radius",RDDx:"Integrated Color Solutions",RDG:"Roland DG",REDM:"REDMS Group",RELI:"Relisys",RGMS:"Rolf Gierling Multitools",RICO:"Ricoh",RNLD:"Edmund Ronald",ROYA:"Royal",RPC:"Ricoh Printing Systems",RTL:"Royal Information Electronics",SAMP:"Sampo",SAMS:"Samsung",SANT:"Jaime Santana Pomares",SCIT:"Scitex",SCRN:"Dainippon Screen",SDP:"Scitex",SEC:"Samsung",SEIK:"Seiko Instruments",SEIk:"Seikosha",SGUY:"ScanGuy.com",SHAR:"Sharp Laboratories",SICC:"International Color Consortium",SONY:"Sony",SPCL:"SpectraCal",STAR:"Star",STC:"Sampo Technology",Scit:"Scitex",Sdp:"Scitex",Sony:"Sony",TALO:"Talon Technology",TAND:"Tandy",TATU:"Tatung",TAXA:"TAXAN America",TDS:"Tokyo Denshi Sekei",TECO:"TECO Information Systems",TEGR:"Tegra",TEKT:"Tektronix",TI:"Texas Instruments",TMKR:"TypeMaker",TOSB:"Toshiba",TOSH:"Toshiba",TOTK:"TOTOKU ELECTRIC",TRIU:"Triumph",TSBT:"Toshiba",TTX:"TTX Computer Products",TVM:"TVM Professional Monitor",TW:"TW Casper",ULSX:"Ulead Systems",UNIS:"Unisys",UTZF:"Utz Fehlau & Sohn",VARI:"Varityper",VIEW:"Viewsonic",VISL:"Visual communication",VIVO:"Vivo Mobile Communication",WANG:"Wang",WLBR:"Wilbur Imaging",WTG2:"Ware To Go",WYSE:"WYSE Technology",XERX:"Xerox",XRIT:"X-Rite",ZRAN:"Zoran",Zebr:"Zebra Technologies",appl:"Apple Computer",bICC:"basICColor",berg:"bergdesign",ceyd:"Integrated Color Solutions",clsp:"MacDermid ColorSpan",ds:"Dainippon Screen",dupn:"DuPont",ffei:"FujiFilm Electronic Imaging",flux:"FluxData",iris:"Scitex",kart:"Scitex",lcms:"Little CMS",lino:"Linotronic",none:"none",ob4d:"Erdt Systems",obic:"Medigraph",quby:"Qubyx Sarl",scit:"Scitex",scrn:"Dainippon Screen",sdp:"Scitex",siwi:"SIWI GRAFIKA",yxym:"YxyMaster"},Qt={scnr:"Scanner",mntr:"Monitor",prtr:"Printer",link:"Device Link",abst:"Abstract",spac:"Color Space Conversion Profile",nmcl:"Named Color",cenc:"ColorEncodingSpace profile",mid:"MultiplexIdentification profile",mlnk:"MultiplexLink profile",mvis:"MultiplexVisualization profile",nkpf:"Nikon Input Device Profile (NON-STANDARD!)"};L($,"icc",[[4,Be],[12,Qt],[40,Object.assign({},Be,Qt)],[48,Be],[80,Be],[64,{0:"Perceptual",1:"Relative Colorimetric",2:"Saturation",3:"Absolute Colorimetric"}],["tech",{amd:"Active Matrix Display",crt:"Cathode Ray Tube Display",kpcd:"Photo CD",pmd:"Passive Matrix Display",dcam:"Digital Camera",dcpj:"Digital Cinema Projector",dmpc:"Digital Motion Picture Camera",dsub:"Dye Sublimation Printer",epho:"Electrophotographic Printer",esta:"Electrostatic Printer",flex:"Flexography",fprn:"Film Writer",fscn:"Film Scanner",grav:"Gravure",ijet:"Ink Jet Printer",imgs:"Photo Image Setter",mpfr:"Motion Picture Film Recorder",mpfs:"Motion Picture Film Scanner",offs:"Offset Lithography",pjtv:"Projection Television",rpho:"Photographic Paper Printer",rscn:"Reflective Scanner",silk:"Silkscreen",twax:"Thermal Wax Printer",vidc:"Video Camera",vidm:"Video Monitor"}]]);class Ve extends Y{static canHandle(e,t,r){return e.getUint8(t+1)===237&&e.getString(t+4,9)==="Photoshop"&&this.containsIptc8bim(e,t,r)!==void 0}static headerLength(e,t,r){let i,s=this.containsIptc8bim(e,t,r);if(s!==void 0)return i=e.getUint8(t+s+7),i%2!=0&&(i+=1),i===0&&(i=4),s+8+i}static containsIptc8bim(e,t,r){for(let i=0;i<r;i++)if(this.isIptcSegmentHead(e,t+i))return i}static isIptcSegmentHead(e,t){return e.getUint8(t)===56&&e.getUint32(t)===943868237&&e.getUint16(t+4)===1028}parse(){let{raw:e}=this,t=this.chunk.byteLength-1,r=!1;for(let i=0;i<t;i++)if(this.chunk.getUint8(i)===28&&this.chunk.getUint8(i+1)===2){r=!0;let s=this.chunk.getUint16(i+3),a=this.chunk.getUint8(i+2),o=this.chunk.getLatin1String(i+5,s);e.set(a,this.pluralizeValue(e.get(a),o)),i+=4+s}else if(r)break;return this.translate(),this.output}pluralizeValue(e,t){return e!==void 0?e instanceof Array?(e.push(t),e):[e,t]:t}}m(Ve,"type","iptc"),m(Ve,"translateValues",!1),m(Ve,"reviveValues",!1),F.set("iptc",Ve),L(N,"iptc",[[0,"ApplicationRecordVersion"],[3,"ObjectTypeReference"],[4,"ObjectAttributeReference"],[5,"ObjectName"],[7,"EditStatus"],[8,"EditorialUpdate"],[10,"Urgency"],[12,"SubjectReference"],[15,"Category"],[20,"SupplementalCategories"],[22,"FixtureIdentifier"],[25,"Keywords"],[26,"ContentLocationCode"],[27,"ContentLocationName"],[30,"ReleaseDate"],[35,"ReleaseTime"],[37,"ExpirationDate"],[38,"ExpirationTime"],[40,"SpecialInstructions"],[42,"ActionAdvised"],[45,"ReferenceService"],[47,"ReferenceDate"],[50,"ReferenceNumber"],[55,"DateCreated"],[60,"TimeCreated"],[62,"DigitalCreationDate"],[63,"DigitalCreationTime"],[65,"OriginatingProgram"],[70,"ProgramVersion"],[75,"ObjectCycle"],[80,"Byline"],[85,"BylineTitle"],[90,"City"],[92,"Sublocation"],[95,"State"],[100,"CountryCode"],[101,"Country"],[103,"OriginalTransmissionReference"],[105,"Headline"],[110,"Credit"],[115,"Source"],[116,"CopyrightNotice"],[118,"Contact"],[120,"Caption"],[121,"LocalCaption"],[122,"Writer"],[125,"RasterizedCaption"],[130,"ImageType"],[131,"ImageOrientation"],[135,"LanguageIdentifier"],[150,"AudioType"],[151,"AudioSamplingRate"],[152,"AudioSamplingResolution"],[153,"AudioDuration"],[154,"AudioOutcue"],[184,"JobID"],[185,"MasterDocumentID"],[186,"ShortDocumentID"],[187,"UniqueDocumentID"],[188,"OwnerID"],[200,"ObjectPreviewFileFormat"],[201,"ObjectPreviewFileVersion"],[202,"ObjectPreviewData"],[221,"Prefs"],[225,"ClassifyState"],[228,"SimilarityIndex"],[230,"DocumentNotes"],[231,"DocumentHistory"],[232,"ExifCameraInfo"],[255,"CatalogSets"]]),L($,"iptc",[[10,{0:"0 (reserved)",1:"1 (most urgent)",2:"2",3:"3",4:"4",5:"5 (normal urgency)",6:"6",7:"7",8:"8 (least urgent)",9:"9 (user-defined priority)"}],[75,{a:"Morning",b:"Both Morning and Evening",p:"Evening"}],[131,{L:"Landscape",P:"Portrait",S:"Square"}]]);function un(n,e,t){const r=t.type===1||t.type===2||t.type===7?1:t.type===3?2:t.type===4||t.type===9||t.type===11||t.type===13?4:t.type===5||t.type===10||t.type===12?8:0;if(r===0)throw new Error(`Unsupported TIFF field type: ${t.type}`);const s=t.count*r<=4?t.valueFieldOffset:n.getUint32(t.valueFieldOffset,e),a=[];for(let o=0;o<t.count;o++){const l=s+o*r;if(l<0||l+r>n.byteLength)throw new Error("Invalid TIFF field offset");if(t.type===1)a.push(n.getUint8(l));else if(t.type===2||t.type===7)a.push(n.getUint8(l));else if(t.type===3)a.push(n.getUint16(l,e));else if(t.type===5){const h=n.getUint32(l+4,e);a.push(h?n.getUint32(l,e)/h:0)}else if(t.type===13)a.push(n.getUint32(l,e));else if(t.type===9)a.push(n.getInt32(l,e));else if(t.type===10){const h=n.getInt32(l+4,e);a.push(h?n.getInt32(l,e)/h:0)}else t.type===11?a.push(n.getFloat32(l,e)):t.type===12?a.push(n.getFloat64(l,e)):a.push(n.getUint32(l,e))}return a}function Zt(n,e,t){if(!Number.isFinite(t)||t<=0||t+2>n.byteLength)throw new Error("Invalid TIFF IFD offset");const r=n.getUint16(t,e);if(t+2+r*12+4>n.byteLength)throw new Error("Corrupt TIFF IFD");const s=new Map;for(let a=0;a<r;a++){const o=t+2+a*12,l=n.getUint16(o,e),h=n.getUint16(o+2,e),c=n.getUint32(o+4,e);s.set(l,{type:h,count:c,valueFieldOffset:o+8})}return s}function H(n,e,t,r){const i=t.get(r);return i?un(n,e,i):[]}const wr=(...n)=>{for(const e of n){if(Array.isArray(e)||ArrayBuffer.isView(e)){const r=Array.from(e).map(Number).filter(i=>Number.isFinite(i)&&i>0);if(r.length)return Math.max(...r);continue}const t=Number(e);if(Number.isFinite(t)&&t>0)return t}return null},fn=n=>{const t=(Array.isArray(n)||ArrayBuffer.isView(n)?Array.from(n).map(Number):[Number(n)]).map(r=>Number.isFinite(r)?r:0);return t.length>=3?[t[0]||0,t[1]||0,t[1]||0,t[2]||0]:t.length===1?[t[0]||0,t[0]||0,t[0]||0,t[0]||0]:[0,0,0,0]};function gn(n,e={}){var k,W,q,X,x,B;if(n.byteLength<8)return null;const t=new DataView(n.buffer,n.byteOffset,n.byteLength),r=t.getUint16(0,!1),i=r===18761;if(!i&&r!==19789||t.getUint16(2,i)!==42)return null;const s=t.getUint32(4,i);let a;try{a=Zt(t,i,s)}catch{return null}const o=new Set([s]);for(const f of H(t,i,a,330))f>0&&o.add(f);let l=null;for(const f of o)try{const C=f===s?a:Zt(t,i,f),D=H(t,i,C,256)[0]||0,S=H(t,i,C,257)[0]||0,w=H(t,i,C,258),R=H(t,i,C,277)[0]||w.length||1,J=w.length===1?new Array(R).fill(w[0]):w,ue=H(t,i,C,259)[0]||1,re=H(t,i,C,262)[0]||0,Q=H(t,i,C,284)[0]||1,V=H(t,i,C,339),rt=V.length===0?new Array(R).fill(1):V.length===1?new Array(R).fill(V[0]):V,ve=H(t,i,C,273),Me=H(t,i,C,279),nt=H(t,i,C,278)[0]||S,it=rt.slice(0,R).every(se=>se===0||se===1);if(!(re===34892&&ue===1&&Q===1&&D>0&&S>0&&R>=3&&J.length>=R&&J.slice(0,R).every(se=>se===16)&&it&&ve.length>0&&ve.length===Me.length))continue;(!l||D*S>l.width*l.height)&&(l={entries:C,width:D,height:S,samplesPerPixel:R,rowsPerStrip:nt,stripOffsets:ve,stripByteCounts:Me,bitsPerSample:J,photometric:re})}catch{continue}if(!l)return null;const{entries:h,width:c,height:d,samplesPerPixel:u,rowsPerStrip:p,stripOffsets:_,stripByteCounts:g,bitsPerSample:y,photometric:b}=l,v=c*d;if(!Number.isSafeInteger(v)||v<=0)return null;const A=new Uint16Array(v*3);for(let f=0;f<_.length;f++){const C=f*p;if(C>=d)break;const D=Math.min(p,d-C),S=D*c*u*2,w=_[f],R=g[f];if(w<0||R<S||w+S>n.byteLength)throw new Error("Invalid LinearRaw DNG strip bounds");const J=C*c*3;if(u===3&&i&&!(n.byteOffset+w&1)){A.set(new Uint16Array(n.buffer,n.byteOffset+w,D*c*3),J);continue}const ue=new DataView(n.buffer,n.byteOffset+w,S);let re=0,Q=J;for(let V=0;V<D*c;V++)A[Q++]=ue.getUint16(re,i),A[Q++]=ue.getUint16(re+2,i),A[Q++]=ue.getUint16(re+4,i),re+=u*2}const T=H(t,i,h,50714),P=((W=(k=e==null?void 0:e.color_data)==null?void 0:k.dng_levels)==null?void 0:W.dng_cblack)||((q=e==null?void 0:e.color_data)==null?void 0:q.cblack_rawpy_style)||(e==null?void 0:e.black_level_per_channel)||(e==null?void 0:e.cblack),O=fn(T.length?T:P),E=H(t,i,h,50717),U=wr(E,(x=(X=e==null?void 0:e.color_data)==null?void 0:X.dng_levels)==null?void 0:x.dng_whitelevel,(B=e==null?void 0:e.color_data)==null?void 0:B.maximum,e==null?void 0:e.white_level)||65535;return{data:A,width:c,height:d,bayerPattern:"",blackLevels:O,whiteLevel:U,metadata:{...e,format:"DNG_LINEAR_RAW_RGB",description:"Uncompressed DNG LinearRaw RGB",linearRawDngDecoder:!0,bitsPerSample:y.slice(0,u),samplesPerPixel:u,photometric:b},isThreePlane:!0,threePlaneTransfer:"linear"}}async function mn(n){const e=gn(n);if(!e)return null;let t={};const r=await bt(),i=new r;try{await i.open(n,{}),t=await i.metadata(!0)}catch(s){console.warn("LinearRaw DNG metadata enrichment failed",s)}finally{try{i.delete?i.delete():i.close()}catch{}}try{const s=await yr.parse(n.buffer);s&&(t={...t,...s})}catch(s){console.warn("exifr parsing failed for LinearRaw DNG",s)}return e.metadata={...t,...e.metadata},e}const _n=async n=>{var i,s,a,o,l,h,c,d,u,p,_;const e=await mn(n);if(e)return e;const t=await bt(),r=new t;try{if(await r.open(n,{}),typeof r.getRawImage!="function")throw new Error("WASM mismatch");let g={};try{g=await r.metadata(!0)}catch(x){console.warn("Metadata error before raw extraction",x)}const y=r.getRawImage(),b=y.data instanceof Uint16Array?y.data:new Uint16Array(y.data);let v={...g};try{const x=await yr.parse(n.buffer);x&&(v={...v,...x})}catch(x){console.warn("exifr parsing failed for RAW buffer",x)}const A=g.filters||((i=g.idata)==null?void 0:i.filters)||0,T=g.colors||((s=g.idata)==null?void 0:s.colors)||0,P=A===0&&T===3,O=A===9;let E=[0,0,0,0],U=null,k;if(r.getBlackLevels)try{const x=r.getBlackLevels();U=x,k=Ft(x)||Ft((a=g.color_data)==null?void 0:a.black_level_model)||void 0;const B=me(k==null?void 0:k.siteBaseLevels);B&&(E=B)}catch(x){console.warn("getBlackLevels binding failed",x)}if(!k){const x=Number((U==null?void 0:U.black)??((o=g.color_data)==null?void 0:o.black)??0)||0,B=me((U==null?void 0:U.cblack)||((l=g.color_data)==null?void 0:l.cblack_rawpy_style)||g.black_level_per_channel||g.cblack||((h=g.color)==null?void 0:h.cblack));if(B){const f=String(y.bayerPattern||g.cfa_pattern||"RGGB").toUpperCase();let C=0;E=[0,1,2,3].map(D=>{const S=f[D],w=S==="R"?0:S==="B"?2:++C===1?1:3;return Math.max(0,x+B[w])})}else E=[x,x,x,x]}const W=Number.isFinite(Number(y.bits))&&Number(y.bits)>0?Math.pow(2,Number(y.bits))-1:null,q=wr(g.white_level,(d=(c=g.color_data)==null?void 0:c.dng_levels)==null?void 0:d.dng_whitelevel,(u=g.color_data)==null?void 0:u.maximum,U==null?void 0:U.maximum,(p=g.color_data)==null?void 0:p.fmaximum,(_=g.color_data)==null?void 0:_.data_maximum,W)||16383,X={data:b,width:y.width,height:y.height,bayerPattern:y.bayerPattern||"",blackLevels:E,blackLevelModel:k,whiteLevel:q,metadata:v,isThreePlane:P,threePlaneTransfer:P?"linear":void 0,isXTrans:O};return zr(X,Vr(X,v)),X}finally{r.delete?r.delete():r.close()}};async function yn(n){if(Lr(n)){const r=await Er(n);if(!r)throw new Error("Sony cRAW HQ decoder did not return image data.");return r.rawImageData}const t=new Uint8Array(n);return _n(t)}function wn(n,e,t,r,i,s=!1,a=!1,o={rows:!1,rowGroup:0,cols:!1,colGroup:0},l=!0,h=1,c=!1){const d=Math.max(0,Math.floor(r.x)),u=Math.max(0,Math.floor(r.y)),p=Math.min(e,Math.ceil(r.x+r.w)),_=Math.min(t,Math.ceil(r.y+r.h)),g=p-d,y=_-u;if(g<=1||y<=1)return null;const b=c&&n.length>=e*t*3;if(b){if(i<0||i>2)return null}else if(i<0||i>3)return null;let v=0,A=0;if(b)v=p-d,A=_-u;else{for(let f=d;f<p;f++)f%2===(i&1)&&v++;for(let f=u;f<_;f++)f%2===i>>1&&A++}if(v<4||A<4)return null;const T=Math.max(1,h),P=Math.floor(2*v/(T+1)),O=Math.floor(2*A/(T+1));if(P<4||O<4)return null;const E=new He(P),U=new He(O),k=new Float32Array(Math.floor(P/2)+1),W=new Float32Array(Math.floor(O/2)+1);let q=0,X=0;const x=new Float32Array(v),B=new Float32Array(A);for(let f=u;f<_;f++){const C=f%2;if(!b&&i>>1!==C||o.rows&&(b?f&1:f>>1&1)!==o.rowGroup)continue;const D=f*e;let S=0;for(let w=d;w<p;w++){if(!b&&w%2!==(i&1)||S>=v)continue;const R=b?n[(D+w)*3+i]:n[D+w],J=Number.isFinite(R)?R:0;x[S++]=J}for(let w=0;w<T;w++){const R=Math.floor(w*P/2);if(R+P>v)break;E.calculateSpectrumWindow(x,R,k,a),q++}if(!l)break}for(let f=d;f<p;f++){const C=f%2;if(!b&&(i&1)!==C||o.cols&&(b?f&1:f>>1&1)!==o.colGroup)continue;let D=0;for(let S=u;S<_;S++){if(!b&&S%2!==i>>1||D>=A)continue;const w=b?n[(S*e+f)*3+i]:n[S*e+f],R=Number.isFinite(w)?w:0;B[D++]=R}for(let S=0;S<T;S++){const w=Math.floor(S*O/2);if(w+O>A)break;U.calculateSpectrumWindow(B,w,W,a),X++}if(!l)break}if(q>0)for(let f=0;f<k.length;f++)k[f]/=q;if(X>0)for(let f=0;f<W.length;f++)W[f]/=X;return{avgHoriz:Array.from(k),avgVert:Array.from(W),sampleCount:{h:q,v:X}}}const er=192*1024*1024,Ye=new Map;let Ge=0,tr=Promise.resolve();const bn=n=>{n.lastAccess=Date.now()},xn=(n=new Set)=>{if(Ge<=er)return;const e=Array.from(Ye.entries()).filter(([t])=>!n.has(t)).sort((t,r)=>t[1].lastAccess-r[1].lastAccess);for(const[t,r]of e){if(Ge<=er)break;Ge-=r.sizeBytes,Ye.delete(t)}},Sn=n=>{var e;return n.data.byteLength+(((e=n.floatData)==null?void 0:e.byteLength)||0)},br=async n=>{const e=Ye.get(n.key);if(e)return bn(e),e.raw;if(!n.file)throw new Error(`Worker cache miss for ${n.key}`);const t=await n.file.arrayBuffer(),r=await yn(t),i={raw:r,sizeBytes:Sn(r),lastAccess:Date.now()};return Ye.set(n.key,i),Ge+=i.sizeBytes,xn(new Set([n.key])),r},kn=async n=>{for(const e of n.files)await br(e);self.postMessage({id:n.id,type:"result",success:!0})},Cn=async n=>{const e=new Map;for(const i of n.files){const s=await br(i);for(const a of n.channelIndices){const o=wn(s.data,s.width,s.height,n.roi,a,!1,n.useHanning,n.interlace,n.doAverage,n.nSegments,!!s.isThreePlane);if(!o)continue;let l=e.get(a);if(!l)l={sumH:new Float32Array(o.avgHoriz.length),sumV:new Float32Array(o.avgVert.length),sampleCount:o.sampleCount,contributingFrames:0},e.set(a,l);else if(l.sumH.length!==o.avgHoriz.length||l.sumV.length!==o.avgVert.length)throw new Error("Stack FFT requires matching frame dimensions and channel layout.");for(let h=0;h<o.avgHoriz.length;h++)l.sumH[h]+=o.avgHoriz[h];for(let h=0;h<o.avgVert.length;h++)l.sumV[h]+=o.avgVert[h];l.sampleCount=o.sampleCount,l.contributingFrames+=1}self.postMessage({id:n.id,type:"progress",completed:1})}const t=Array.from(e.entries()).map(([i,s])=>({ch:i,sumH:s.sumH,sumV:s.sumV,sampleCount:s.sampleCount,contributingFrames:s.contributingFrames})),r=[];for(const i of t)r.push(i.sumH.buffer,i.sumV.buffer);self.postMessage({id:n.id,type:"result",success:!0,processedFrames:n.files.length,results:t},r)};self.onmessage=n=>{const e=n.data;tr=tr.then(async()=>{try{e.type==="preload"?await kn(e):await Cn(e)}catch(t){self.postMessage({id:e.id,type:"result",success:!1,error:(t==null?void 0:t.message)||"Power spectrum worker failed"})}})};
