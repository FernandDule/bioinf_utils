# Bioinf Utils

## Описание
Набор утилит для обработки биоинформатических данных (FASTQ, FASTA, BLAST).

## Структура пакета
- `main.py` – функции-утилиты.
- `helpers/fastq.py` – вспомогательные функции для фильтрации FASTQ.
- `helpers/bio_files.py` – вспомогательные функции для чтения/записи FASTA и парсинга BLAST.

## Доступные утилиты

### 1. `filter_fastq`
Фильтрует FASTQ-файл по GC-составу, длине ридов и среднему качеству. Отфильтрованные риды сохраняются в папку `./filtered`.
* `input_fastq` (str): Путь к входному FASTQ-файлу.
* `output_fastq` (str): Имя выходного файла.
* `gc_bounds` (tuple | int | float, по умолч. `(0, 100)`): Границы GC-состава в %. Если передано одно число, оно считается верхней границей.
* `length_bounds` (tuple | int, по умолч. `(0, 2**32)`): Границы длины рида. Если передано одно число, оно считается верхней границей.
* `quality_threshold` (int | float, по умолч. `0`): Минимальный средний балл качества (шкала Phred33).

### 2. `convert_multiline_fasta_to_oneline`
Преобразует многострочный FASTA-файл в однострочный формат.
* `input_fasta` (str): Путь к входному FASTA-файлу.
* `output_fasta` (str | None, по умолч. `None`): Путь к выходному файлу. Если не указан, файл создается в той же директории с суффиксом `_oneline.fasta`.

### 3. `parse_blast_output`
Извлекает названия best hits из текстового отчёта BLAST и сохраняет их, отсортированными по алфавиту.
* `input_file` (str): Путь к текстовому файлу с результатами BLAST.
* `output_file` (str | None, по умолч. `None`): Путь к выходному файлу. Если не указан, файл создается в той же директории с суффиксом `_selected_blast_queries.txt`.

## Вызов (примеры использования)

```python
from main import filter_fastq, convert_multiline_fasta_to_oneline, parse_blast_output

# 1. Фильтрация FASTQ
filter_fastq(
    input_fastq="data/input.fastq",
    output_fastq="cleaned.fastq",
    gc_bounds=(40, 60),
    length_bounds=50,
    quality_threshold=20
)

# 2. Преобразование FASTA
convert_multiline_fasta_to_oneline("data/multiline.fasta")

# 3. Парсинг BLAST
parse_blast_output("data/blast_results.txt")