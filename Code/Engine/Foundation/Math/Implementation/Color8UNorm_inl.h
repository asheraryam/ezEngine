#pragma once

EZ_FORCE_INLINE ezColorLinearUB::ezColorLinearUB(ezUInt8 R, ezUInt8 G, ezUInt8 B, ezUInt8 A /* = 255*/) : ezColorUnsignedByteBase(R, G, B, A)
{
}

inline ezColorLinearUB::ezColorLinearUB(const ezColor& color)
{
  *this = color;
}

inline void ezColorLinearUB::operator=(const ezColor& color)
{
  r = ezMath::ColorFloatToByte(color.r);
  g = ezMath::ColorFloatToByte(color.g);
  b = ezMath::ColorFloatToByte(color.b);
  a = ezMath::ColorFloatToByte(color.a);
}

inline ezColor ezColorLinearUB::ToLinearFloat() const
{
  const float f = 1.0f / 255.0f;

  return ezColor(r * f, g * f, b * f, a * f);
}

// *****************

EZ_FORCE_INLINE ezColorGammaUB::ezColorGammaUB(ezUInt8 R, ezUInt8 G, ezUInt8 B, ezUInt8 A) : ezColorUnsignedByteBase(R, G, B, A)
{
}

inline ezColorGammaUB::ezColorGammaUB(const ezColor& color)
{
  *this = color;
}

inline void ezColorGammaUB::operator=(const ezColor& color)
{
  const ezVec3 gamma = ezColor::LinearToGamma(ezVec3(color.r, color.g, color.b));

  r = ezMath::ColorFloatToByte(gamma.x);
  g = ezMath::ColorFloatToByte(gamma.y);
  b = ezMath::ColorFloatToByte(gamma.z);
  a = ezMath::ColorFloatToByte(color.a);
}

inline ezColor ezColorGammaUB::ToLinearFloat() const
{
  ezVec3 gamma;
  gamma.x = ezMath::ColorByteToFloat(r);
  gamma.y = ezMath::ColorByteToFloat(g);
  gamma.z = ezMath::ColorByteToFloat(b);

  const ezVec3 linear = ezColor::GammaToLinear(gamma);

  return ezColor(linear.x, linear.y, linear.z, a * (1.0f / 255.0f));
}


