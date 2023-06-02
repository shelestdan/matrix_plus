#include "s21_matrix_oop.h"

S21Matrix::S21Matrix() : rows_(0), cols_(0), matrix_(nullptr) {}

S21Matrix::~S21Matrix() { DestroyMatrix(); }

S21Matrix::S21Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
  AllocateMemory();
}

S21Matrix::S21Matrix(const S21Matrix &other)
    : rows_(other.rows_), cols_(other.cols_) {
  AllocateMemory();
  CopyMatrix(other);
}

S21Matrix::S21Matrix(S21Matrix &&other)
    : rows_(other.rows_), cols_(other.cols_), matrix_(other.matrix_) {
  other.rows_ = 0;
  other.cols_ = 0;
  other.matrix_ = nullptr;
}

S21Matrix &S21Matrix::operator=(const S21Matrix &other) {
  if (this != &other) {
    DestroyMatrix();
    rows_ = other.rows_;
    cols_ = other.cols_;
    AllocateMemory();
    CopyMatrix(other);
  }
  return *this;
}

S21Matrix &S21Matrix::operator+=(const S21Matrix &other) {
  SumMatrix(other);
  return *this;
}

S21Matrix &S21Matrix::operator-=(const S21Matrix &other) {
  SubMatrix(other);
  return *this;
}

S21Matrix &S21Matrix::operator*=(const S21Matrix &other) {
  MulMatrix(other);
  return *this;
}

S21Matrix &S21Matrix::operator*=(const double num) {
  MulNumber(num);
  return *this;
}

S21Matrix S21Matrix::operator+(const S21Matrix &other) {
  S21Matrix result(*this);
  result += other;
  return result;
}

S21Matrix S21Matrix::operator-(const S21Matrix &other) {
  S21Matrix result(*this);
  result -= other;
  return result;
}

S21Matrix S21Matrix::operator*(const S21Matrix &other) {
  S21Matrix result(*this);
  result *= other;
  return result;
}

S21Matrix S21Matrix::operator*(const double num) {
  S21Matrix result(*this);
  result *= num;
  return result;
}

double &S21Matrix::operator()(int row, int col) const {
  return matrix_[row][col];
}

void S21Matrix::SetRowsAndCols(int rows, int cols) {
  if (rows <= 0 || cols <= 0) {
    throw std::out_of_range("New matrix size can't be less then 1");
  }

  double **new_mtx = new double *[rows];
  for (int row = 0; row < rows; row++) {
    new_mtx[row] = new double[cols];

    for (int col = 0; col < cols; col++) {
      if (row < rows_ && col < cols_) {
        new_mtx[row][col] = matrix_[row][col];
      } else {
        new_mtx[row][col] = 0.0;
      }
    }
  }

  DestroyMatrix();

  rows_ = rows;
  cols_ = cols;

  matrix_ = new_mtx;
}

void S21Matrix::SumMatrix(const S21Matrix &other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::out_of_range("Matrix dimensions do not match!");
  }

  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      matrix_[i][j] += other.matrix_[i][j];
    }
  }
}

void S21Matrix::SubMatrix(const S21Matrix &other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::out_of_range("Error: Matrices have different sizes.");
  }

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] -= other.matrix_[i][j];
    }
  }
}

void S21Matrix::MulNumber(const double num) {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] *= num;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix &other) {
  if (cols_ != other.rows_) {
    throw std::out_of_range(
        "Error: Matrices cannot be multiplied due to incompatible sizes.");
  }

  S21Matrix result(rows_, other.cols_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < other.cols_; j++) {
      double sum = 0;
      for (int k = 0; k < cols_; k++) {
        sum += matrix_[i][k] * other.matrix_[k][j];
      }
      result.matrix_[i][j] = sum;
    }
  }

  DestroyMatrix();
  CopyMatrix(result);
  result.DestroyMatrix();
}

S21Matrix operator*(const double num, const S21Matrix &other) {
  S21Matrix result(other.rows_, other.cols_);
  for (int i = 0; i < other.rows_; i++) {
    for (int j = 0; j < other.cols_; j++) {
      result.matrix_[i][j] = other.matrix_[i][j] * num;
    }
  }
  return result;
}

bool S21Matrix::EqMatrix(const S21Matrix &other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    return false;
  }

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      if (!AreEqual(matrix_[i][j], other.matrix_[i][j], E)) {
        return false;
      }
    }
  }

  return true;
}

bool S21Matrix::AreEqual(double a, double b, double epsilon) {
  return std::abs(a - b) <= epsilon;
}

void S21Matrix::DestroyMatrix() {
  for (int i = 0; i < rows_; i++) {
    delete[] matrix_[i];
  }
  delete[] matrix_;
  matrix_ = nullptr;
  rows_ = 0;
  cols_ = 0;
}

void S21Matrix::AllocateMemory() {
  matrix_ = new double *[rows_];
  for (int i = 0; i < rows_; i++) {
    double *rowPtr = new double[cols_];
    std::fill(rowPtr, rowPtr + cols_, 0.0);
    matrix_[i] = rowPtr;
  }
}

void S21Matrix::CopyMatrix(const S21Matrix &other) {
  DestroyMatrix();
  rows_ = other.rows_;
  cols_ = other.cols_;
  matrix_ = new double *[rows_];
  for (int row = 0; row < rows_; row++) {
    matrix_[row] = new double[cols_];
    for (int col = 0; col < cols_; col++) {
      matrix_[row][col] = other.matrix_[row][col];
    }
  }
}

S21Matrix S21Matrix::MinorMatrix(int i_row, int j_col) const {
  if (rows_ < 2 || cols_ < 2) {
    throw std::invalid_argument(
        "Cannot calculate minor matrix for a matrix with "
        "dimensions less than 2x2");
  }

  S21Matrix minor(rows_ - 1, cols_ - 1);

  for (int i = 0; i < rows_ - 1; ++i) {
    for (int j = 0; j < cols_ - 1; ++j) {
      minor(i, j) = (*this)(i < i_row ? i : i + 1, j < j_col ? j : j + 1);
    }
  }

  return minor;
}

S21Matrix S21Matrix::Transpose() const {
  S21Matrix TransposedMatrix(cols_, rows_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      TransposedMatrix(j, i) = matrix_[i][j];
    }
  }
  return TransposedMatrix;
}

S21Matrix S21Matrix::CalcComplements() {
  S21Matrix com(rows_, cols_);

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      S21Matrix minor = MinorMatrix(i, j);
      double det = minor.Determinant();
      double sign = ((i + j) % 2 == 0) ? 1.0 : -1.0;
      com(i, j) = sign * det;
    }
  }

  return com;
}

S21Matrix S21Matrix::InverseMatrix() {
  double det = Determinant();
  if (std::abs(det) < E) {
    throw std::invalid_argument("Matrix is singular.");
  }

  S21Matrix com = CalcComplements();
  S21Matrix result = com.Transpose();
  result.MulNumber(1 / det);

  return result;
}

bool S21Matrix::operator==(const S21Matrix &other) { return EqMatrix(other); }

double S21Matrix::Determinant() const {
  if (rows_ != cols_) {
    throw std::invalid_argument(
        "Error: Determinant is only defined for square matrices.");
  }

  if (rows_ == 1) {
    return matrix_[0][0];
  }

  double det = 0;
  for (int j = 0; j < cols_; j++) {
    S21Matrix minor = MinorMatrix(0, j);
    double sign = (j % 2 == 0) ? 1.0 : -1.0;
    det += sign * matrix_[0][j] * minor.Determinant();
  }

  return det;
}
