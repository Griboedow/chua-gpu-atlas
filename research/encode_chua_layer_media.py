"""Encode the published layer cards as lossless WebP and an optional seekable MP4."""
import argparse
import hashlib
import json
import shutil
import subprocess
from pathlib import Path

from PIL import Image, ImageChops, features

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/chua-layer-cards'
ASSETS = ROOT / 'docs/assets'


def source_frames():
    data = json.loads((OUT / 'cards.json').read_text(encoding='utf-8'))
    paths = [ASSETS / f'chua-layer-cards/{id}.png' for id in data['frame_order']]
    frames = []
    for path in paths:
        with Image.open(path) as image:
            frames.append(image.convert('RGB'))
    assert len(frames) == len(data['frame_durations_ms']) == data['frame_count']
    assert all(list(frame.size) == data['image_size'] for frame in frames)
    return data, paths, frames


def verify_webp(path, frames, durations):
    with Image.open(path) as encoded:
        assert encoded.n_frames == len(frames)
        assert encoded.info['loop'] == 0
        elapsed = 0
        for index, (source, duration) in enumerate(zip(frames, durations)):
            encoded.seek(index)
            encoded.load()
            assert encoded.size == source.size
            assert encoded.info['duration'] == duration
            assert encoded.info['timestamp'] == elapsed
            assert ImageChops.difference(source, encoded.convert('RGB')).getbbox() is None, index
            elapsed += duration
    return dict(decoded_frames=len(frames), duration_ms=elapsed, all_rgb_pixels_equal_to_png=True)


def encode_webp(force=False):
    data, paths, frames = source_frames()
    target = ASSETS / 'chua-layer-cards.webp'
    current = False
    if target.exists() and not force:
        try:
            verify_webp(target, frames, data['frame_durations_ms'])
            current = True
        except (AssertionError, OSError, KeyError):
            pass
    if not current:
        print('Encoding lossless WebP', flush=True)
        frames[0].save(target, save_all=True, append_images=frames[1:],
                       duration=data['frame_durations_ms'], loop=0,
                       lossless=True, quality=80, method=4, minimize_size=True)
    checks = verify_webp(target, frames, data['frame_durations_ms'])
    gif = ASSETS / 'chua-layer-cards.gif'
    report = dict(path=target.relative_to(ROOT).as_posix(), bytes=target.stat().st_size,
                  gif_bytes=gif.stat().st_size, lossless=True, quality=80, method=4,
                  minimize_size=True, libwebp_version=features.version('webp'), **checks)
    print('WebP:', report['bytes'], 'bytes; all source pixels and durations preserved', flush=True)
    return report, data, paths, frames


def encode_video(ffmpeg, data, frames):
    width, height = data['image_size']
    durations = data['frame_durations_ms']
    assert width % 2 == height % 2 == 0
    assert all(duration > 0 and duration % 1000 == 0 for duration in durations)
    output = ASSETS / 'chua-layer-cards.mp4'
    scratch = ROOT / '.publication/layer-media'
    scratch.mkdir(parents=True, exist_ok=True)
    command = [str(ffmpeg), '-y', '-hide_banner', '-loglevel', 'warning',
               '-f', 'rawvideo', '-pixel_format', 'rgb24', '-video_size', f'{width}x{height}',
               '-framerate', '1', '-i', 'pipe:0', '-an', '-c:v', 'libx264',
               '-preset', 'slow', '-tune', 'stillimage', '-crf', '23', '-pix_fmt', 'yuv420p',
               '-vf', 'scale=1280:-2:flags=lanczos',
               '-r', '1', '-g', '8', '-movflags', '+faststart', str(output)]
    with (scratch / 'video-encode.log').open('wb') as log:
        with subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.DEVNULL, stderr=log) as process:
            for id, frame, duration in zip(data['frame_order'], frames, durations):
                pixels = frame.tobytes()
                for _ in range(duration // 1000):
                    process.stdin.write(pixels)
                print('Video:', id, flush=True)
            process.stdin.close()
            assert process.wait() == 0, 'FFmpeg failed; see .publication/layer-media/video-encode.log'
    # Decode the whole stream, not just its header. The output is lossy H.264.
    subprocess.run([str(ffmpeg), '-v', 'error', '-i', str(output), '-f', 'null', '-'], check=True)
    return dict(path=output.relative_to(ROOT).as_posix(), bytes=output.stat().st_size,
                codec='H.264', pixel_format='yuv420p', crf=23, preset='slow', tune='stillimage',
                source_dimensions=[width, height], encoded_dimensions=[1280, 2*round(height*1280/width/2)],
                frames_per_second=1, frames=sum(durations)//1000, duration_ms=sum(durations),
                keyframe_interval=8, faststart=True, full_decode_passed=True,
                ffmpeg_version=subprocess.check_output([str(ffmpeg), '-version'], text=True).splitlines()[0])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--video', action='store_true', help='Also create an H.264 video with native player seeking.')
    parser.add_argument('--ffmpeg', type=Path, help='FFmpeg executable; otherwise use PATH.')
    parser.add_argument('--force', action='store_true', help='Re-encode WebP even if its decoded frames already match.')
    args = parser.parse_args()
    webp, data, paths, frames = encode_webp(args.force)
    report = dict(webp=webp, source_frame_order=data['frame_order'],
                  source_png_sha256={path.relative_to(ROOT).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},
                  source_durations_ms=data['frame_durations_ms'])
    if args.video:
        ffmpeg = args.ffmpeg or shutil.which('ffmpeg')
        if not ffmpeg:
            raise SystemExit('For video, install FFmpeg or pass --ffmpeg PATH.')
        report['mp4'] = encode_video(ffmpeg, data, frames)
    (OUT / 'media.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8', newline='\n')


if __name__ == '__main__':
    main()
