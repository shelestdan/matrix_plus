#ifndef SRC_S21_MATRIX_OOP_H_
#define SRC_S21_MATRIX_OOP_H_

#include <cmath>
#include <iostream>
#define E 1e-7

class S21Matrix {
 public:
  S21Matrix();
  ~S21Matrix();
  S21Matrix(int rows, int cols);
  S21Matrix(const S21Matrix &other);
  S21Matrix(S21Matrix &&other);

  bool EqMatrix(const S21Matrix &other);
  void SumMatrix(const S21Matrix &other);
  void SubMatrix(const S21Matrix &other);
  void MulNumber(const double num);
  void MulMatrix(const S21Matrix &other);

  double &operator()(int row, int col) const;
  friend S21Matrix operator*(const double num, const S21Matrix &other);
  bool operator==(const S21Matrix &other);
  S21Matrix operator+(const S21Matrix &other);
  S21Matrix operator-(const S21Matrix &other);
  S21Matrix operator*(const S21Matrix &other);
  S21Matrix operator*(const double num);
  S21Matrix &operator=(const S21Matrix &other);
  S21Matrix &operator+=(const S21Matrix &other);
  S21Matrix &operator-=(const S21Matrix &other);
  S21Matrix &operator*=(const S21Matrix &other);
  S21Matrix &operator*=(const double num);

  S21Matrix Transpose() const;
  S21Matrix CalcComplements();
  S21Matrix InverseMatrix();
  double Determinant() const;
  void SetRows(int rows) { SetRowsAndCols(rows, cols_); }
  void SetColumns(int cols) { SetRowsAndCols(rows_, cols); }
  int GetRows() { return rows_; }
  int GetColumns() { return cols_; }
  void SetRowsAndCols(int rows, int cols);

 private:
  int rows_, cols_;
  double **matrix_;
  void DestroyMatrix();
  void AllocateMemory();
  void CopyMatrix(const S21Matrix &other);
  bool AreEqual(double a, double b, double epsilon = E);

  S21Matrix MinorMatrix(int i_rows, int j_cols) const;
};

#endif  // SRC_S21_MATRIX_OOP_H_