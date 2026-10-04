// The bundled libraylib.a was built with an msvcrt-based MinGW, where stat() resolves to
// stat64i32. UCRT-based MinGW toolchains only export _stat64i32, so the missing name is
// provided here. CMake only adds this file when it is needed.
extern "C" {
struct _stat64i32;
int _stat64i32(const char* path, struct _stat64i32* buffer);
int stat64i32(const char* path, struct _stat64i32* buffer) { return _stat64i32(path, buffer); }
}
