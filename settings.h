#pragma once

// ===================== SETTINGS =====================

// true  -> приближённая формула
// false -> интеграл (как в оригинале)
constexpr bool USE_APPROX = false;

// защита от NaN / inf / слишком больших значений
constexpr bool USE_CUTS = true;

// минимальное значение (для защиты)
constexpr double MIN_VAL = 1e-30;

// директория с файлами
constexpr const char* DATA_DIR = ".";