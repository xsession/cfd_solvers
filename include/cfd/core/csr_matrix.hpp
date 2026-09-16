#pragma once

#include "cfd/core/parallel.hpp"

#include <algorithm>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace cfd::core {

class CsrMatrix {
public:
    CsrMatrix() = default;
    CsrMatrix(std::size_t rows,
              std::size_t cols,
              std::vector<std::size_t> row_offsets,
              std::vector<std::size_t> column_indices,
              std::vector<double> values)
        : rows_(rows),
          cols_(cols),
          row_offsets_(std::move(row_offsets)),
          column_indices_(std::move(column_indices)),
          values_(std::move(values)) {
        validate();
    }

    [[nodiscard]] std::size_t rows() const noexcept { return rows_; }
    [[nodiscard]] std::size_t cols() const noexcept { return cols_; }
    [[nodiscard]] std::size_t nonzeros() const noexcept { return values_.size(); }
    [[nodiscard]] const std::vector<std::size_t>& row_offsets() const noexcept { return row_offsets_; }
    [[nodiscard]] const std::vector<std::size_t>& column_indices() const noexcept { return column_indices_; }
    [[nodiscard]] const std::vector<double>& values() const noexcept { return values_; }

    void multiply(std::span<const double> x, std::span<double> y) const {
        if (x.size() != cols_ || y.size() != rows_) {
            throw std::invalid_argument("CSR multiply size mismatch");
        }
        parallel_for(rows_, [&](std::size_t row) {
            double sum = 0.0;
            for (std::size_t k = row_offsets_[row]; k < row_offsets_[row + 1U]; ++k) {
                sum += values_[k] * x[column_indices_[k]];
            }
            y[row] = sum;
        });
    }

    [[nodiscard]] std::vector<double> diagonal() const {
        std::vector<double> result(rows_, 0.0);
        parallel_for(rows_, [&](std::size_t row) {
            for (std::size_t k = row_offsets_[row]; k < row_offsets_[row + 1U]; ++k) {
                if (column_indices_[k] == row) {
                    result[row] = values_[k];
                    break;
                }
            }
        });
        return result;
    }

private:
    std::size_t rows_{};
    std::size_t cols_{};
    std::vector<std::size_t> row_offsets_;
    std::vector<std::size_t> column_indices_;
    std::vector<double> values_;

    void validate() const {
        if (row_offsets_.size() != rows_ + 1U) {
            throw std::invalid_argument("CSR row-offset size mismatch");
        }
        if (column_indices_.size() != values_.size()) {
            throw std::invalid_argument("CSR column/value size mismatch");
        }
        if (row_offsets_.empty() || row_offsets_.front() != 0U || row_offsets_.back() != values_.size()) {
            throw std::invalid_argument("CSR row offsets are invalid");
        }
        for (std::size_t row = 0; row < rows_; ++row) {
            if (row_offsets_[row] > row_offsets_[row + 1U]) {
                throw std::invalid_argument("CSR row offsets must be monotone");
            }
            std::size_t previous = 0U;
            bool have_previous = false;
            for (std::size_t k = row_offsets_[row]; k < row_offsets_[row + 1U]; ++k) {
                const std::size_t column = column_indices_[k];
                if (column >= cols_) throw std::invalid_argument("CSR column index out of range");
                if (have_previous && column <= previous) {
                    throw std::invalid_argument("CSR columns must be strictly increasing in each row");
                }
                previous = column;
                have_previous = true;
            }
        }
    }
};

class CsrBuilder {
public:
    CsrBuilder(std::size_t rows, std::size_t cols) : rows_(rows), cols_(cols), entries_(rows) {}

    void add(std::size_t row, std::size_t column, double value) {
        if (row >= rows_ || column >= cols_) throw std::out_of_range("CSR builder index out of range");
        entries_[row].emplace_back(column, value);
    }

    [[nodiscard]] CsrMatrix build(double drop_tolerance = 0.0) const {
        if (drop_tolerance < 0.0) throw std::invalid_argument("CSR drop tolerance must be non-negative");
        std::vector<std::size_t> row_offsets(rows_ + 1U, 0U);
        std::vector<std::size_t> columns;
        std::vector<double> values;

        for (std::size_t row = 0; row < rows_; ++row) {
            auto row_entries = entries_[row];
            std::sort(row_entries.begin(), row_entries.end(), [](const auto& a, const auto& b) {
                return a.first < b.first;
            });
            std::size_t i = 0U;
            while (i < row_entries.size()) {
                const std::size_t column = row_entries[i].first;
                double value = 0.0;
                do {
                    value += row_entries[i].second;
                    ++i;
                } while (i < row_entries.size() && row_entries[i].first == column);
                if (value > drop_tolerance || value < -drop_tolerance) {
                    columns.push_back(column);
                    values.push_back(value);
                }
            }
            row_offsets[row + 1U] = values.size();
        }
        return CsrMatrix(rows_, cols_, std::move(row_offsets), std::move(columns), std::move(values));
    }

private:
    std::size_t rows_{};
    std::size_t cols_{};
    std::vector<std::vector<std::pair<std::size_t, double>>> entries_;
};

} // namespace cfd::core
